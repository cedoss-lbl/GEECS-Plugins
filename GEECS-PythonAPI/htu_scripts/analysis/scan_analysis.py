import os
import time
import inspect
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from progressbar import ProgressBar
from typing import Union, NamedTuple, Any, Optional
from geecs_api.interface import GeecsDatabase, api_error
from geecs_api.devices.geecs_device import GeecsDevice
from geecs_api.devices.HTU.diagnostics.cameras import Camera
from geecs_api.tools.distributions.binning import unsupervised_binning, BinningResults
from geecs_api.tools.images.scan_images import ScanImages
from geecs_api.tools.scans.scan_data import ScanData
from geecs_api.tools.images.filtering import FiltersParameters
from geecs_api.tools.interfaces.exports import load_py, save_py
from geecs_api.tools.interfaces.prompts import text_input
# from geecs_api.tools.distributions.fit_utility import fit_distribution


def scan_analysis(scan_data: ScanData, device: Union[GeecsDevice, str], variable: str, camera: Union[int, Camera, str],
                  blind_loads: bool = False, store_images: bool = True, save: bool = False) \
        -> tuple[Optional[Path], dict[str, Any]]:
    analyses: list[dict[str, Any]] = []

    scan_images = ScanImages(scan_data, camera)
    device_name: str = device.get_name() if isinstance(device, GeecsDevice) else device
    key_data = scan_data.data_dict[device_name]

    # scan parameters & binning
    measured: BinningResults = unsupervised_binning(key_data[variable], key_data['shot #'])

    Expected = NamedTuple('Expected', start=float, end=float, steps=int, shots=int, setpoints=np.ndarray, indexes=list)
    steps: int = 1 + round((float(scan_data.scan_info['End']) - float(scan_data.scan_info['Start']))
                           / float(scan_data.scan_info['Step size']))
    expected = Expected(start=float(scan_data.scan_info['Start']),
                        end=float(scan_data.scan_info['End']),
                        steps=steps,
                        shots=int(scan_data.scan_info['Shots per step']),
                        setpoints=np.linspace(float(scan_data.scan_info['Start']),
                                              float(scan_data.scan_info['End']),
                                              steps),
                        indexes=[np.arange(p * steps, (p+1) * steps) for p in range(steps)])

    matching = \
        all([inds.size == expected.shots for inds in measured.indexes]) and (len(measured.indexes) == expected.steps)
    if not matching:
        api_error.warning(f'Observed data binning does not match expected scan parameters (.ini)',
                          f'Function "{inspect.stack()[0][3]}"')

    # list images for each step
    def build_file_name(shot: int):
        return scan_images.image_folder / \
            f'Scan{scan_data.get_tag().number:03d}_{scan_images.camera_name}_{shot:03d}.png'

    if matching:
        indexes = expected.indexes
        setpoints = expected.setpoints
    else:
        indexes = measured.indexes
        setpoints = measured.avg_x

    paths: list[list[Path]] = \
        [[build_file_name(ind+1) for ind in np.sort(inds) if build_file_name(ind+1).is_file()] for inds in indexes]

    # run image analyses
    analysis_files: list[Path] = []

    with ProgressBar(max_value=len(paths)) as pb:
        for it, (step_paths, step_val) in enumerate(zip(paths, setpoints)):
            # check if analysis exists
            keep: Union[str, bool] = False
            save_dir: Path = scan_data.get_analysis_folder() / f'Step_{it+1}'
            analysis_file: Union[Path, str] = save_dir / 'profiles_analysis.dat'

            analyze: str = 'y'
            if analysis_file.is_file():
                analyze = text_input(f'\nRe-run the analysis (step "{step_val}")? : ',
                                     accepted_answers=['y', 'yes', 'n', 'no'])

            # run/load analysis
            analysis: dict[str, Any] = {}
            if (analyze.lower()[0] == 'y') or (not analysis_file.is_file()):
                print(f'\nAnalyzing step "{step_val}"...')
                if not save_dir.is_dir():
                    os.makedirs(save_dir)

                scan_images.set_save_folder(save_dir)
                analysis_file, analysis = \
                    scan_images.run_analysis_with_checks(images=step_paths,
                                                         initial_filtering=FiltersParameters(com_threshold=0.66),
                                                         plots=True, store_images=store_images, save=save)

            if not analysis:
                print('Loading analysis...')
                analysis, analysis_file = load_py(analysis_file)
                keep = blind_loads
                if not blind_loads:
                    ScanImages.render_image_analysis(analysis['average_analysis'], tag='average_image', block=True)

            if not analysis:
                continue  # skip

            if not analysis_file:
                analysis_file = ''
            analysis_files.append(analysis_file)

            if keep:
                keep = 'y'
            else:
                keep = text_input(f'Add this analysis to the overall screen scan analysis? : ',
                                  accepted_answers=['y', 'yes', 'n', 'no'])
            if keep.lower()[0] == 'n':
                continue

            keep = text_input(f'Add this analysis to the overall screen scan analysis? : ',
                              accepted_answers=['y', 'yes', 'n', 'no'])
            if keep.lower()[0] == 'n':
                continue

            print('Collecting analysis summary...')
            analyses.append(analysis)

            pb.increment()
            time.sleep(0.01)

    # export to .dat
    data_dict: dict[str, Any] = {
        'indexes': indexes,
        'setpoints': setpoints,
        'analysis_files': analysis_files,
        'analyses': analyses,
        'device_name': device_name,
        'scan_folder': scan_images.scan.get_folder(),
        'camera_name': scan_images.camera_name,
        'pos_short_names': pos_short_names,
        'pos_long_names': pos_long_names}
    export_file_path = scan_data.get_analysis_folder() / f'steering_analysis_{device_name}'
    save_py(file_path=export_file_path, data=data_dict, as_bulk=False)
    print(f'Data exported to:\n\t{export_file_path}.dat')

    return export_file_path, data_dict


def render_scan_analysis(data_dict: dict[str, Any], physical_units: bool = True,
                         show_xy: bool = True, show_fwhm: bool = True, show_deltas: bool = True,
                         xy_metric: str = 'mean', fwhm_metric: str = 'mean', deltas_metric: str = 'mean'):
    """
    metric:     'mean', 'median'
    """
    x_axis: np.ndarray = data_dict['setpoints']
    analyses: list[dict[str, Any]] = data_dict['analyses']
    n_rows: int = sum([show_xy, show_fwhm, show_deltas])

    for pos, label in zip(analyses[0]['positions']['short_names'], analyses[0]['positions']['long_names']):
        fig, axs = plt.subplots(ncols=1, nrows=n_rows,
                                figsize=(ScanImages.fig_size[0], ScanImages.fig_size[1] * 1.5),
                                sharex='col', sharey='row')
        # X(var), Y(var)
        axs[0].plot(x_axis, analyses['max_mean_pos_pix'][:, 0], '.k', markersize=10)
        opt = data_dict['beam_analysis']['x_fit']['opt']
        opt_sign = '-' if opt[1] < 0 else ''
        axs[0, it].plot(x_axis, analyses['x_fit']['fit'], 'gray',
                        label=rf"$X \simeq {opt[0]} \cdot dx {opt_sign}{abs(opt[1])} \sigma$")

        axs[0, it].plot(x_axis, analyses['max_mean_pos_pix'][:, 0],
                        '.k', markersize=10)
        opt = data_dict['beam_analysis']['x_fit']['opt']
        opt_sign = '-' if opt[1] < 0 else ''
        axs[0, it].plot(x_axis, analyses['x_fit']['fit'], 'gray',
                        label=rf"$X \simeq {opt[0]} \cdot dx {opt_sign}{abs(opt[1])} \sigma$")

        axs[0, it].legend(loc='best', prop={'size': 8})
        axs[0, it].set_xticks([])
        axs[0, it].set_title(pos_long_names[it])

        # X(dy), Y(dy)
        axs[1, it].fill_between(
            x_axis,
            f_deltas * (beam_analysis[f'{pos}_deltas_means'][:, 0] - beam_analysis[f'{pos}_deltas_stds'][:, 0]),
            f_deltas * (beam_analysis[f'{pos}_deltas_means'][:, 0] + beam_analysis[f'{pos}_deltas_stds'][:, 0]),
            label=r'$D_y \pm \sigma$', color='m', alpha=0.33)
        axs[1, it].plot(x_axis, f_deltas * beam_analysis[f'{pos}_deltas_avg_imgs'][:, 0], 'ob-',
                        label=r'$D_y$ $(\mu_{image})$', linewidth=1, markersize=3)
        axs[1, it].legend(loc='best', prop={'size': 8})
        axs[1, it].set_xticks([])

        # FWHM X
        axs[2, it].fill_between(
            x_axis,
            f_fwhms * beam_analysis[f'{pos}_fwhm_means'][:, 1] - beam_analysis[f'{pos}_fwhm_stds'][:, 1],
            f_fwhms * beam_analysis[f'{pos}_fwhm_means'][:, 1] + beam_analysis[f'{pos}_fwhm_stds'][:, 1],
            label=r'$FWHM_x \pm \sigma$', color='y', alpha=0.33)
        axs[2, it].plot(x_axis, f_fwhms * beam_analysis[f'{pos}_fwhm_means'][:, 1], 'og-',
                        label=r'$FWHM_x$ $(\mu_{image})$', linewidth=1, markersize=3)
        axs[2, it].legend(loc='best', prop={'size': 8})
        axs[2, it].set_xticks([])

        # FWHM Y
        axs[3, it].fill_between(
            x_axis,
            f_fwhms * beam_analysis[f'{pos}_fwhm_means'][:, 0] - beam_analysis[f'{pos}_fwhm_stds'][:, 0],
            f_fwhms * beam_analysis[f'{pos}_fwhm_means'][:, 0] + beam_analysis[f'{pos}_fwhm_stds'][:, 0],
            label=r'$FWHM_y \pm \sigma$ [$\mu$m]', color='y', alpha=0.33)
        axs[3, it].plot(x_axis, f_fwhms * beam_analysis[f'{pos}_fwhm_means'][:, 0], 'og-',
                        label=r'$FWHM_y$ $(\mu_{image})$', linewidth=1, markersize=3)
        axs[3, it].legend(loc='best', prop={'size': 8})
        axs[3, it].set_xlabel('Screen')
        axs[3, it].set_xticks(x_axis, screen_labels)

        axs[0, 0].set_ylabel(f'X-Offsets [{units_deltas}]')
        axs[1, 0].set_ylabel(f'Y-Offsets [{units_deltas}]')
        axs[2, 0].set_ylabel(f'X-FWHM [{units_fwhms}]')
        axs[3, 0].set_ylabel(f'Y-FWHM [{units_fwhms}]')

        # set matching vertical limits for deltas/FWHMs
        y_lim = (min(axs[0, 0].get_ylim()[0], axs[1, 0].get_ylim()[0]),
                 max(axs[0, 0].get_ylim()[1], axs[1, 0].get_ylim()[1]))
        [axs[0, j].set_ylim(y_lim) for j in range(len(pos_short_names))]
        [axs[1, j].set_ylim(y_lim) for j in range(len(pos_short_names))]

        y_lim = (min(axs[2, 0].get_ylim()[0], axs[3, 0].get_ylim()[0]),
                 max(axs[2, 0].get_ylim()[1], axs[3, 0].get_ylim()[1]))
        [axs[2, j].set_ylim(y_lim) for j in range(len(pos_short_names))]
        [axs[3, j].set_ylim(y_lim) for j in range(len(pos_short_names))]

        if save_dir:
            save_path = save_dir / 'beam_analysis.png'
            plt.savefig(save_path, dpi=300)

        plt.show(block=True)


if __name__ == '__main__':
    GeecsDevice.exp_info = GeecsDatabase.collect_exp_info('Undulator')

    # _base = Path(r'C:\Users\GuillaumePlateau\Documents\LBL\Data')
    _base: Path = Path(r'Z:\data')

    _base_tag = (2023, 4, 13, 26)
    _key_device = 'U_S4H'
    _camera_tag = 'A2'

    _scan_data = ScanData(tag=_base_tag, experiment_base_path=_base / 'Undulator')

    # data preview
    # _key_data = _scan_data.data_dict[_key_device]
    # _bins: BinningResults = unsupervised_binning(_key_data['Current'], _key_data['shot #'])
    #
    # plt.figure()
    # for x, _ind in zip(_bins.avg_x, _bins.indexes):
    #     plt.plot(x * np.ones(_ind.shape), _ind, '.', alpha=0.3)
    # plt.xlabel('Current [A]')
    # plt.ylabel('Indexes')
    # plt.show(block=True)

    # run analysis
    _export_file_path, _data_dict = steering_scan_analysis(_scan_data, _key_device, _camera_tag)

    # open analysis
    # _analysis_file = Path(r'Z:\data\Undulator\Y2023\04-Apr\23_0413\analysis\Scan026\steering_analysis_U_S4H.dat')
    # _data_dict, _ = load_py(_analysis_file, as_dict=True)

    print('done')
