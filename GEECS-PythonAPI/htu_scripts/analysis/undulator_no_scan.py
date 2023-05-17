import cv2
import os
import re
import math
import numpy as np
import screeninfo
from pathlib import Path
import scipy.ndimage as simg
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.transform import hough_ellipse
from skimage.feature import canny
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from progressbar import ProgressBar
from typing import Optional, Any, Union
from geecs_api.api_defs import SysPath
from geecs_api.tools.images.batches import list_images
from geecs_api.devices.geecs_device import api_error
import geecs_api.tools.images.ni_vision as ni
from geecs_api.tools.scans.scan import Scan
from geecs_api.devices.geecs_device import GeecsDevice
from geecs_api.devices.HTU.diagnostics.cameras import Camera
from geecs_api.tools.images.filtering import clip_hot_pixels, filter_image
from geecs_api.tools.images.spot import spot_analysis, fwhm


class UndulatorNoScan:
    fig_size = (int(round(screeninfo.get_monitors()[0].width / 540. * 10) / 10),
                int(round(screeninfo.get_monitors()[0].height / 450. * 10) / 10))

    def __init__(self, scan: Scan, camera: Union[int, Camera, str], angle: Optional[int] = None):
        """
        Container for data analysis of a set of images collected at an undulator station.

        scan (Scan object): analysis object for the relevant scan
        camera (int | Camera | str), either of:
            - station number (1-9)
            - Camera object
            - GEECS device name of the relevant camera
            - relevant screen shorthand label (U1 - U9)
        angle (int): rotation angle to apply (multiples of +/-90 deg only). Ignored if camera object is provided.
        """

        self.scan: Scan = scan
        self.camera: Optional[Camera] = None
        if angle is None:
            self.camera_r90: int = 0
        else:
            self.camera_r90: int = int(round(angle / 90.))

        if isinstance(camera, Camera):
            self.camera = camera
            self.camera_name: str = camera.get_name()
            self.camera_roi: Optional[np.ndarray] = camera.roi
            self.camera_r90 = camera.rot_90
        elif isinstance(camera, str) and (camera in Camera.ROIs):
            self.camera_name = camera
            self.camera_roi = np.array(Camera.ROIs[camera])
        elif isinstance(camera, str) and re.match(r'(U[1-9]|A[1-3]|Rad2)', camera):
            self.camera_name = Camera.name_from_label(camera)
            self.camera_roi = np.array(Camera.ROIs[self.camera_name])
        elif isinstance(camera, int) and (1 <= camera <= 9):
            self.camera_name = Camera.name_from_label(f'U{camera}')
            self.camera_roi = np.array(Camera.ROIs[self.camera_name])
        else:
            self.camera_name = camera
            self.camera_roi = None

        self.camera_label: str = Camera.label_from_name(self.camera_name)

        self.image_folder: SysPath = self.scan.get_folder() / self.camera_name
        self.image_analyses: Optional[list[dict[str, Any]]] = None
        self.analyses_summary: Optional[dict[str, Any]] = None

        self.average_image: Optional[np.ndarray] = None
        self.average_analysis: Optional[dict[str, Any]] = None

    @staticmethod
    def find_roi(image: np.ndarray, threshold: Optional[float] = None, plots: bool = False):
        roi = None
        roi_box = np.array([0, image.shape[1] - 1, 0, image.shape[0]])

        try:
            # filter and smooth
            blur = clip_hot_pixels(image, median_filter_size=2, threshold_factor=3)
            blur = simg.gaussian_filter(blur, sigma=5.)

            # threshold
            if threshold is None:
                counts, bins = np.histogram(blur, bins=10)
                threshold = bins[np.where(counts == np.max(counts))[0][0] + 1]
            bw = closing(blur > threshold, square(3))

            # remove artifacts connected to image border
            cleared = clear_border(bw)
            # cleared = bw

            # label image regions
            label_image = label(cleared)
            areas = [box.area for box in regionprops(label_image)]
            roi = regionprops(label_image)[areas.index(max(areas))]
            roi_box = roi.bbox

        except Exception:
            pass

        finally:
            if plots:
                fig, ax = plt.subplots(figsize=UndulatorNoScan.fig_size)
                ax.imshow(image)
                rect = mpatches.Rectangle((roi_box[1], roi_box[0]), roi_box[3] - roi_box[1], roi_box[2] - roi_box[0],
                                          fill=False, edgecolor='red', linewidth=2)
                ax.add_patch(rect)
                ax.set_axis_off()
                plt.tight_layout()
                plt.show(block=True)

        return roi_box, roi

    def read_image_as_float(self, image_path: SysPath) -> np.ndarray:
        image = ni.read_imaq_image(image_path)
        if isinstance(self.camera_roi, np.ndarray) and (self.camera_roi.size >= 4):
            image = image[self.camera_roi[-2]:self.camera_roi[-1], self.camera_roi[0]:self.camera_roi[1]]
        image = np.rot90(image, self.camera_r90)
        return image.astype('float64')

    def analyze_images(self, n_images: int = 0, contrast: float = 2, hp_median: int = 2, hp_threshold: float = 3.,
                       denoise_cycles: int = 0, gauss_filter: float = 5., com_threshold: float = 0.5,
                       plots: bool = False, save_dir: Optional[Path] = None, skip_ellipse: bool = True):
        paths = list_images(self.image_folder, n_images, '.png')
        # paths = paths[24:27]  # tmp
        paths = paths[:3]  # tmp
        analyses: list[dict[str, Any]] = []

        # run analysis
        if paths:
            try:
                with ProgressBar(max_value=len(paths)) as pb:
                    # analyze one at a time until valid image found to initialize average
                    for skipped, image_path in enumerate(paths):
                        analysis = self._analyze_image(image_path, contrast, hp_median, hp_threshold, denoise_cycles,
                                                       gauss_filter, com_threshold, skip_ellipse)
                        analysis['image_path'] = image_path

                        if analysis['is_valid']:
                            analysis = UndulatorNoScan.profiles_analysis(analysis)
                            if plots:
                                self.plot_save_image_analysis(analysis, block=True)
                            avg_image: np.ndarray = analysis['image_raw'].copy()
                            break

                    analyses.append(analysis)
                    pb.increment()

                    # analyze the rest of the images
                    if len(paths) > (skipped + 1):
                        for it, image_path in enumerate(paths[skipped+1:]):
                            analysis = self._analyze_image(image_path, contrast, hp_median, hp_threshold,
                                                           denoise_cycles, gauss_filter, com_threshold, skip_ellipse)
                            analysis['image_path'] = image_path

                            if analysis['is_valid']:
                                # profiles
                                analysis = UndulatorNoScan.profiles_analysis(analysis)
                                if plots:
                                    self.plot_save_image_analysis(analysis, block=True)

                                # average
                                alpha = 1.0 / (it + 2)
                                beta = 1.0 - alpha
                                avg_image = cv2.addWeighted(analysis['image_raw'], alpha, avg_image, beta, 0.0)

                            # collect
                            analyses.append(analysis)
                            pb.increment()

                # analyze the average image
                self.average_analysis = self._analyze_image(avg_image, contrast, hp_median, hp_threshold,
                                                            denoise_cycles, gauss_filter, com_threshold, skip_ellipse)
                if self.average_analysis['is_valid']:
                    self.average_analysis = \
                        UndulatorNoScan.profiles_analysis(self.average_analysis)
                    if plots:
                        self.plot_save_image_analysis(self.average_analysis, block=True)

                self.image_analyses = analyses
                self.average_image = avg_image
                self.analyses_summary = self._summarize_image_analyses()
                self.analyses_summary['targets'] = self._target_analysis()

            except Exception as ex:
                api_error.error(str(ex), 'Failed to analyze image')
                pass

    def _analyze_image(self, image: Union[SysPath, np.ndarray], contrast: float = 2,
                       hp_median: int = 2, hp_threshold: float = 3., denoise_cycles: int = 0,
                       gauss_filter: float = 5., com_threshold: float = 0.5,
                       skip_ellipse: bool = True) -> dict[str, Any]:
        # raw mage
        if isinstance(image, np.ndarray):
            image_raw = image
        else:
            image_raw = self.read_image_as_float(image)

        # initial filtering
        analysis = filter_image(image_raw, hp_median, hp_threshold, denoise_cycles, gauss_filter, com_threshold)

        # check image (contrast)
        counts, bins = np.histogram(analysis['image_denoised'],
                                    bins=max(10, 2 * int(np.std(analysis['image_blurred']))))
        analysis['bkg_level']: float = bins[np.where(counts == np.max(counts))[0][0]]
        analysis['is_valid']: bool = np.std(bins) > contrast * analysis['bkg_level']

        # stop if low-contrast image
        if not analysis['is_valid']:
            return analysis

        try:
            # beam edges
            image_edges: np.ndarray = canny(analysis['image_thresholded'], sigma=3)
            image_edges = closing(image_edges.astype(float), square(3))
            image_edges = clear_border(image_edges)

            # boxed image
            bw = closing(analysis['image_blurred'] > 0.5 * np.max(analysis['image_blurred']), square(3))
            cleared = clear_border(bw)
            label_image = label(cleared)
            areas = [box.area for box in regionprops(label_image)]
            roi = regionprops(label_image)[areas.index(max(areas))]

            roi = np.array([roi.bbox[1], roi.bbox[3], roi.bbox[0], roi.bbox[2]])  # left, right, top, bottom
            width_gain = int(round((roi[1] - roi[0]) * 0.1))
            width_gain = min(min(roi[0], width_gain), min(image_edges.shape[1] - 1 - roi[1], width_gain))
            height_gain = int(round((roi[3] - roi[2]) * 0.1))
            height_gain = min(min(roi[2], height_gain), min(image_edges.shape[0] - 1 - roi[3], height_gain))
            analysis['roi'] = np.array([roi[0] - width_gain, roi[1] + width_gain,
                                        roi[2] - height_gain, roi[3] + height_gain])

            pos_box = np.array([(analysis['roi'][2] + analysis['roi'][3]) / 2.,
                                (analysis['roi'][0] + analysis['roi'][1]) / 2.])
            pos_box = np.round(pos_box).astype(int)
            analysis['position_box'] = tuple(pos_box)  # i, j

            # update edges image
            image_edges[:, :analysis['roi'][0]] = 0
            image_edges[:, analysis['roi'][1]+1:] = 0
            image_edges[:analysis['roi'][2], :] = 0
            image_edges[analysis['roi'][3]+1:, :] = 0
            analysis['image_edges'] = image_edges

            # ellipse fit (min_size: min of major axis, max_size: max of minor axis)
            if not skip_ellipse:
                h_ellipse = hough_ellipse(image_edges, min_size=10, max_size=60, accuracy=8)
                h_ellipse.sort(order='accumulator')
                best = list(h_ellipse[-1])
                yc, xc, a, b = (int(round(x)) for x in best[1:5])
                ellipse = (yc, xc, a, b, best[5])
                analysis['ellipse'] = ellipse  # (i, j, major, minor, orientation)

                # cy, cx = ellipse_perimeter(*ellipse)
                pos_ellipse = np.array([yc, xc])
                analysis['position_ellipse'] = tuple(pos_ellipse)  # i, j

        except Exception:
            pass

        return analysis

    def _target_analysis(self) -> dict[str, Any]:
        positions = ['max', 'com', 'box', 'ellipse']
        target_analysis: dict[str, Any] = {}

        try:
            target_analysis['target_um_pix'] = \
                float(GeecsDevice.exp_info['devices'][self.camera_name]['SpatialCalibration']['defaultvalue'])
            target_analysis['target_xy'] = \
                (float(GeecsDevice.exp_info['devices'][self.camera_name]['Target.Y']['defaultvalue']),
                 float(GeecsDevice.exp_info['devices'][self.camera_name]['Target.X']['defaultvalue']))
        except Exception:
            target_analysis['target_um_pix'] = 1.
            target_analysis['target_xy'] = (0., 0.)

        if self.average_analysis and self.analyses_summary:
            for pos in positions:
                if f'position_{pos}' in self.average_analysis:
                    target_analysis[f'avg_img_{pos}_delta'] = \
                        (np.array(self.average_analysis[f'position_{pos}']) - np.array(target_analysis['target_xy'])) \
                        * target_analysis['target_um_pix'] / 1000.

                if f'scan_pos_{pos}' in self.analyses_summary:
                    target_analysis[f'scan_pos_{pos}_delta'] = \
                        (np.array(self.analyses_summary[f'scan_pos_{pos}']) - np.array(target_analysis['target_xy'])) \
                        * target_analysis['target_um_pix'] / 1000.
                    target_analysis[f'target_deltas_{pos}_mean'] = \
                        np.array([np.mean(target_analysis[f'scan_pos_{pos}_delta'], axis=0)])
                    target_analysis[f'target_deltas_{pos}_std'] = \
                        np.array([np.std(target_analysis[f'scan_pos_{pos}_delta'], axis=0)])

        return target_analysis

    def _summarize_image_analyses(self) -> dict[str, Any]:
        positions = ['max', 'com', 'box', 'ellipse']
        summary: dict[str, Any] = {}

        for pos in positions:
            summary[f'scan_pos_{pos}'] = \
                np.array([analysis[f'position_{pos}'] for analysis in self.image_analyses
                          if analysis is not None and f'position_{pos}' in analysis])
            summary[f'scan_pos_{pos}_fwhm_x'] = \
                np.array([fwhm(analysis[f'opt_x_{pos}'][3]) for analysis in self.image_analyses
                          if analysis is not None and f'opt_x_{pos}' in analysis])
            summary[f'scan_pos_{pos}_fwhm_y'] = \
                np.array([fwhm(analysis[f'opt_y_{pos}'][3]) for analysis in self.image_analyses
                          if analysis is not None and f'opt_y_{pos}' in analysis])

            if summary[f'scan_pos_{pos}'].any() and \
                    summary[f'scan_pos_{pos}_fwhm_x'].any() and \
                    summary[f'scan_pos_{pos}_fwhm_y'].any():
                summary[f'mean_pos_{pos}'] = np.mean(summary[f'scan_pos_{pos}'], axis=0)
                summary[f'mean_pos_{pos}_fwhm_x'] = np.mean(summary[f'scan_pos_{pos}_fwhm_x'])
                summary[f'mean_pos_{pos}_fwhm_y'] = np.mean(summary[f'scan_pos_{pos}_fwhm_y'])

                summary[f'std_pos_{pos}'] = np.std(summary[f'scan_pos_{pos}'], axis=0)
                summary[f'std_pos_{pos}_fwhm_x'] = np.std(summary[f'scan_pos_{pos}_fwhm_x'])
                summary[f'std_pos_{pos}_fwhm_y'] = np.std(summary[f'scan_pos_{pos}_fwhm_y'])

        return summary

    @staticmethod
    def profiles_analysis(analysis: dict[str, Any]) -> dict[str, Any]:
        try:
            pos_max: tuple[int, int] = analysis['position_max']
            pos_com: tuple[int, int] = analysis['position_com']
            pos_box: tuple[int, int] = analysis['position_box']
            if 'position_ellipse' in analysis:
                pos_ell: tuple[int, int] = analysis['position_ellipse']
            else:
                pos_ell = (0, 0)

            positions = [(pos_max[0], pos_max[1], 'max'),
                         (pos_com[0], pos_com[1], 'com'),
                         (pos_box[0], pos_box[1], 'box')]
            if 'position_ellipse' in analysis:
                positions.append((pos_ell[0], pos_ell[1], 'ell'))

            labels = ['maximum', 'center of mass', 'box', 'ellipse']

            profiles = spot_analysis(analysis['image_raw'], positions,
                                     x_window=(analysis['roi'][0], analysis['roi'][1]),
                                     y_window=(analysis['roi'][2], analysis['roi'][3]))

            for k, v in profiles.items():
                analysis[k] = v

            analysis['positions'] = positions
            analysis['positions_labels'] = labels

        except Exception as ex:
            print(ex)
            pass

        return analysis

    def plot_save_image_analysis(self, analysis: dict[str, Any], block: bool = False):
        try:
            image_path: Path = Path(analysis['image_path'])
            camera_folder: str = image_path.parts[-2]
            shot_name: str = image_path.name.split(".")[0].split("_")[-1]
            save_folder: Path = Path(self.scan.get_analysis_folder()) / camera_folder
            if not save_folder.is_dir():
                os.makedirs(save_folder)

            fig = plt.figure(figsize=(UndulatorNoScan.fig_size[0] * 1.5, UndulatorNoScan.fig_size[1]))
            grid = plt.GridSpec(2, 4, hspace=0.3, wspace=0.3)
            ax_i = fig.add_subplot(grid[:, :2])
            ax_x = fig.add_subplot(grid[0, 2:])
            ax_y = fig.add_subplot(grid[1, 2:], sharey=ax_x)

            # raw image
            edges = np.where(analysis['image_edges'] != 0)
            ax_i.imshow(analysis['image_raw'], cmap='hot', aspect='equal', origin='upper')
            ax_i.scatter(edges[1], edges[0], s=0.3, c='b', alpha=0.3)

            # roi box
            roi_color = '--y'
            roi_line = 0.66
            ax_i.plot(analysis['roi'][0] * np.ones((analysis['roi'][3] - analysis['roi'][2] + 1)),
                      np.arange(analysis['roi'][2], analysis['roi'][3] + 1), roi_color, linewidth=roi_line)  # left
            ax_i.plot(analysis['roi'][1] * np.ones((analysis['roi'][3] - analysis['roi'][2] + 1)),
                      np.arange(analysis['roi'][2], analysis['roi'][3] + 1), roi_color, linewidth=roi_line)  # right
            ax_i.plot(np.arange(analysis['roi'][0], analysis['roi'][1] + 1),
                      analysis['roi'][2] * np.ones((analysis['roi'][1] - analysis['roi'][0] + 1)),
                      roi_color, linewidth=roi_line)  # top
            ax_i.plot(np.arange(analysis['roi'][0], analysis['roi'][1] + 1),
                      analysis['roi'][3] * np.ones((analysis['roi'][1] - analysis['roi'][0] + 1)),
                      roi_color, linewidth=roi_line)  # bottom

            # ellipse
            if 'ellipse' in analysis:
                ell = mpatches.Ellipse(xy=(analysis['ellipse'][1], analysis['ellipse'][0]),
                                       width=2*analysis['ellipse'][3], height=2*analysis['ellipse'][2],
                                       angle=math.degrees(analysis['ellipse'][4]))
                ax_i.add_artist(ell)
                ell.set_alpha(0.66)
                ell.set_edgecolor('g')
                ell.set_linewidth(1)
                # noinspection PyArgumentList
                ell.set_fill(False)
                # ax_i.plot(analysis['position_ellipse'][1], analysis['position_ellipse'][0], 'k.', markersize=2)

            # max
            ax_i.axvline(analysis['position_max'][1], color='gray', linestyle='--', linewidth=0.5)
            ax_i.axhline(analysis['position_max'][0], color='gray', linestyle='--', linewidth=0.5)
            ax_i.plot(analysis['position_max'][1], analysis['position_max'][0], 'k.', markersize=3)
            if ax_i.xaxis_inverted():
                ax_i.invert_xaxis()
            if not ax_i.yaxis_inverted():
                ax_i.invert_yaxis()

            # lineouts
            ax_x.plot(analysis['axis_x_max'], analysis['data_x_max'], 'b-', label='data')
            ax_x.plot(analysis['axis_x_max'], analysis['fit_x_max'], 'm-',
                      label=f'FWHM: {fwhm(analysis["opt_x_max"][3]):.1f}')
            ax_x.legend(loc='best', prop={'size': 8})

            ax_y.plot(analysis['axis_y_max'], analysis['data_y_max'], 'b-', label='data')
            ax_y.plot(analysis['axis_y_max'], analysis['fit_y_max'], 'm-',
                      label=f'FWHM: {fwhm(analysis["opt_y_max"][3]):.1f}')
            ax_y.legend(loc='best', prop={'size': 8})

            # plot & save
            file_name: str = f'max_profiles_{shot_name}.png'
            image_path: Path = save_folder / file_name
            plt.savefig(image_path, dpi=300)
            plt.show(block=False)

            # ___________________________________________
            if 'positions' in analysis:
                profiles_fig_size = (UndulatorNoScan.fig_size[0] * 1.5,
                                     UndulatorNoScan.fig_size[1] * math.ceil(len(analysis['positions']) / 3))
                _, axs = plt.subplots(ncols=3, nrows=len(analysis['positions']), figsize=profiles_fig_size,
                                      sharex='col', sharey='col')
                for it, pos in enumerate(analysis['positions']):
                    axs[it, 0].imshow(analysis['image_raw'])
                    axs[it, 0].axvline(pos[1], color='r', linestyle='--', linewidth=1)
                    axs[it, 0].axhline(pos[0], color='r', linestyle='--', linewidth=1)
                    axs[it, 0].plot(pos[1], pos[0], '.k', markersize=2)
                    axs[it, 0].set_ylabel(analysis['positions_labels'][it])
                    axs[it, 1].plot(analysis[f'axis_x_{pos[2]}'], analysis[f'data_x_{pos[2]}'], 'b-', label='data')
                    axs[it, 1].plot(analysis[f'axis_x_{pos[2]}'], analysis[f'fit_x_{pos[2]}'], 'm-', label='fit(x)')
                    axs[it, 1].legend(loc='best', prop={'size': 8})
                    axs[it, 2].plot(analysis[f'axis_y_{pos[2]}'], analysis[f'data_y_{pos[2]}'], 'b-', label='data')
                    axs[it, 2].plot(analysis[f'axis_y_{pos[2]}'], analysis[f'fit_y_{pos[2]}'], 'm-', label='fit(y)')
                    axs[it, 2].legend(loc='best', prop={'size': 8})

                file_name: str = f'all_profiles_{shot_name}.png'
                image_path: Path = save_folder / file_name
                plt.savefig(image_path, dpi=300)
                plt.show(block=block)

        except Exception as ex:
            print(ex)
            pass


if __name__ == '__main__':
    _base = Path(r'C:\Users\GuillaumePlateau\Documents\LBL\Data')
    # _base = Path(r'Z:\data')
    _folder = _base / r'Undulator\Y2023\05-May\23_0509\scans\Scan028'

    _scan = Scan(_folder, ignore_experiment_name=True)
    images = UndulatorNoScan(_scan, 'U7', angle=0)
    images.analyze_images(contrast=1, hp_median=2, hp_threshold=3., denoise_cycles=0,
                          gauss_filter=5., com_threshold=0.66, plots=True, skip_ellipse=True)
    print('done')