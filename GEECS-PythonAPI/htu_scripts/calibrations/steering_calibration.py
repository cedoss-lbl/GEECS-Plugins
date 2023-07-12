import time
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from typing import Optional, Any
from geecs_api.api_defs import SysPath
from geecs_api.interface import GeecsDatabase
from geecs_api.devices.geecs_device import GeecsDevice
from geecs_api.devices.HTU.transport import Steering
from geecs_api.devices.HTU.diagnostics.screens import Screen
from geecs_api.devices.HTU.diagnostics import EBeamDiagnostics
from geecs_api.devices.HTU.diagnostics.cameras import Camera
from geecs_api.tools.distributions.binning import unsupervised_binning
from geecs_api.tools.interfaces.prompts import text_input
from geecs_api.tools.images.scan_images import ScanImages
from geecs_api.tools.scans.scan_data import ScanData


def steering_calibration(steering_magnets: list[int], screens: list[str],
                         currents: list[tuple[float, float, float, int]],
                         backgrounds: int = 0, live_analysis: bool = True) -> list[dict[str, Any]]:
    """
    Performs current sweeps for each steering magnets and both planes (horizontal/vertical),
    while monitoring e-beam on calibration screens.

    screens:    tuple of screens' shorthand labels e.g., ("P1", "A1")
    currents:   tuple of tuples, where each array consists of
                    "start value" (float),
                    "end value" (float),
                    "step size" (float),
                    "shots per step" (int)
    """
    if not steering_magnets or not screens or not currents:
        return []

    # steering magnets
    magnets = [Steering(n) for n in steering_magnets]
    for magnet in magnets:
        magnet.subscribe_var_values()

    # observation screens
    e_diagnostics = EBeamDiagnostics()
    imagers = [e_diagnostics.imagers[s] for s in screens]

    def close_all():
        for _magnet in magnets:
            _magnet.close()
        e_diagnostics.close()

    # sweeps
    # ------------------------------------------------------------------------------
    sweep_scans: list[dict[str, Any]] = []

    success, scans = sweep_magnet(magnets, currents, imagers[0].camera,
                                  screen_out=None, screen_in=imagers[0].screen, backgrounds=backgrounds)
    sweep_scans += scans
    if not success:
        close_all()
        return sweep_scans

    if len(steering_magnets) > 1:
        for previous_screen, next_screen in zip(imagers[:-1], imagers[1:]):
            success, scans = sweep_magnet(magnets, currents, next_screen.camera,
                                          screen_out=previous_screen.screen, screen_in=next_screen.screen,
                                          backgrounds=backgrounds)
            sweep_scans += scans
            if not success:
                close_all()
                return sweep_scans

    # clear screens
    set_screens(screen_out=imagers[-1].screen, screen_in=None)

    time.sleep(1.)
    close_all()

    return sweep_scans


def sweep_magnet(magnets: list[Steering], currents: list[tuple[float, float, float, int]],
                 camera: Camera, screen_out: Optional[Screen], screen_in: Optional[Screen], backgrounds: int = 0)\
        -> tuple[bool, list[dict[str, Any]]]:
    success: bool = True

    if not set_screens(screen_out, screen_in):
        return False, []

    # background
    if backgrounds > 0:
        print(f'Collecting {backgrounds} background images...')
        camera.save_local_background(n_images=backgrounds)

    scans: list[dict[str, Any]] = []
    for magnet, current in zip(magnets, currents):
        # -------------------
        plane = 'horizontal'
        scan_folder, scan_number, success = \
            run_scan(magnet, plane, current, Camera.label_from_name(camera.get_name()), 300.)
        if success:
            scans.append({'folder': scan_folder,
                          'number': scan_number,
                          'magnet': magnet,
                          'plane': plane,
                          'camera': camera})
        else:
            return False, scans

        # -------------------
        plane = 'vertical'
        scan_folder, scan_number, success = \
            run_scan(magnet, plane, current, Camera.label_from_name(camera.get_name()), 300.)
        if success:
            scans.append({'folder': scan_folder,
                          'number': scan_number,
                          'magnet': magnet,
                          'plane': plane,
                          'camera': camera})
        else:
            return False, scans

    return success, scans


def run_scan(magnet: Steering, plane: str, setpoints: tuple[float, float, float, int],
             screen_tag: str, timeout: float = 300.) -> tuple[Path, int, bool]:
    success: bool = True

    print(f'Starting "S{magnet.get_name()[-1]}" {plane} scan on "{screen_tag}"...')
    while True:
        next_folder, next_scan, accepted, timed_out = magnet.scan_current(plane, *setpoints, timeout=timeout)
        if accepted and not timed_out:
            break
        elif not accepted:
            repeat = text_input(f'Command not accepted. Try again? : ',
                                accepted_answers=['y', 'yes', 'n', 'no'])
            if repeat.lower()[0] == 'n':
                success = False
                break
        else:
            repeat = text_input(f'Scan timed out. Run again (y), ignore (i), stop (s)? : ',
                                accepted_answers=['y', 'yes', 'i', 'ignore', 's', 'stop'])
            if repeat.lower()[0] == 'y':
                continue

            if repeat.lower()[0] == 'i':
                break

            if repeat.lower()[0] == 's':
                success = False
                break

    return next_folder, next_scan, success


def set_screens(screen_out: Optional[Screen], screen_in: Optional[Screen]) -> bool:
    success: bool = True

    # remove screen
    if screen_out is not None:
        print(f'Removing {screen_out.var_alias} ({screen_out.controller.get_name()})...')
        while True:
            screen_out.remove()
            if not screen_out.is_inserted():
                break
            else:
                print(f'Failed to remove {screen_out.var_alias} ({screen_out.controller.get_name()})')
                repeat = text_input(f'Try again? : ', accepted_answers=['y', 'yes', 'n', 'no'])
                if repeat.lower()[0] == 'n':
                    success = False
                    break

    # check
    if not success:
        return success

    # insert screen
    if screen_in is not None:
        print(f'Inserting {screen_in.var_alias} ({screen_in.controller.get_name()})...')
        while True:
            screen_in.insert()
            if screen_in.is_inserted():
                break
            else:
                print(f'Failed to insert {screen_in.var_alias} ({screen_in.controller.get_name()})')
                repeat = text_input(f'Try again? : ', accepted_answers=['y', 'yes', 'n', 'no'])
                if repeat.lower()[0] == 'n':
                    success = False
                    break

    return success


if __name__ == '__main__':
    # base_path = Path(r'C:\Users\GuillaumePlateau\Documents\LBL\Data')
    base_path: Path = Path(r'Z:\data')

    is_local = (str(base_path)[0] == 'C')
    if not is_local:
        GeecsDevice.exp_info = GeecsDatabase.collect_exp_info('Undulator')

    _sweep_scans = steering_calibration(steering_magnets=[1, 2],
                                        screens=['DP', 'P1'],
                                        currents=[(-1., 1., 2., 3),
                                                  (-1., 1., 2., 3)],
                                        backgrounds=10)

    print('done')
