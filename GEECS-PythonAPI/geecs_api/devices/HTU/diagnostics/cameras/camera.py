from __future__ import annotations
import os
import glob
import time
import shutil
import cv2
from typing import Any
from datetime import datetime as dtime
from geecs_api.api_defs import VarAlias, SysPath
from geecs_api.devices.geecs_device import GeecsDevice


class Camera(GeecsDevice):
    def __init__(self, device_name: str):
        super().__init__(device_name)

        # self.gui_path: SysPath = GeecsDevice.exp_info['GUIs'][device_name]

        self.__variables = {VarAlias('BackgroundPath'): (None, None),
                            VarAlias('localsavingpath'): (None, None),
                            VarAlias('exposure'): (None, None),
                            VarAlias('triggerdelay'): (None, None)}
        self.build_var_dicts(tuple(self.__variables.keys()))
        self.var_bkg_path: str = self.var_names_by_index.get(0)[0]
        self.var_save_path: str = self.var_names_by_index.get(1)[0]
        self.var_exposure: str = self.var_names_by_index.get(2)[0]

        # self.register_cmd_executed_handler()
        # self.register_var_listener_handler()

    def get_variables(self):
        return self.__variables

    def interpret_value(self, var_alias: VarAlias, val_string: str) -> Any:
        if var_alias in self.__variables.keys():
            if (var_alias == VarAlias('BackgroundPath')) or (var_alias == VarAlias('localsavingpath')):
                return val_string
            else:
                try:
                    return float(val_string)
                except Exception:
                    return 0.
        else:
            return val_string

    def save_local_background(self, n_images: int = 15):
        # background folders
        avg_bkg_folder: SysPath = os.path.join(self.data_root_path, 'backgrounds', f'{self.get_name()}')
        if not os.path.isdir(avg_bkg_folder):
            os.makedirs(avg_bkg_folder)

        bkg_imgs_folder: SysPath = os.path.join(avg_bkg_folder, 'tmp_images')
        if os.path.isdir(bkg_imgs_folder):
            shutil.rmtree(bkg_imgs_folder, ignore_errors=True)
        os.makedirs(bkg_imgs_folder)
        self.set(self.var_save_path, value=bkg_imgs_folder, exec_timeout=10., sync=True)

        # file name
        stamp = dtime.now()
        file_name = stamp.strftime('%y%m%d_%H%M') + f'_{self.get_name()}_{n_images}_images_avg_bkg.png'

        # save images
        if n_images > 0:
            self.set('save', 'on', exec_timeout=5.)
            time.sleep(n_images + 2)
            self.set('save', 'off', exec_timeout=5.)

        # average image
        Camera.calculate_average_image(self, bkg_imgs_folder, avg_bkg_folder, file_name, n_images)

    def save_background(self, exec_timeout: float = 30.0):
        # background folder
        next_scan_folder, _ = self.next_scan_folder()
        # bkg_target_folder: SysPath = os.path.join(next_scan_folder, f'{self.get_name()}_Background')

        # exp_bkg_folder: SysPath = os.path.join(self.data_root_path, 'backgrounds')
        # if not os.path.isdir(exp_bkg_folder):
        #     os.makedirs(exp_bkg_folder)

        bkg_target_folder: SysPath = os.path.join(self.data_root_path, 'backgrounds', f'{self.get_name()}')
        if not os.path.isdir(bkg_target_folder):
            os.makedirs(bkg_target_folder)

        stamp = dtime.now()
        file_name = stamp.strftime('%y%m%d_%H%M') + f'_{self.get_name()}_avg_bkg.png'

        # save images
        GeecsDevice.run_no_scan(monitoring_device=self, comment=f'{self.get_name()}: background collection',
                                timeout=exec_timeout)

        # average image
        bkg_source_path: SysPath = os.path.join(next_scan_folder, self.get_name())
        Camera.calculate_average_image(self, bkg_source_path, bkg_target_folder, file_name)

    @staticmethod
    def save_multiple_backgrounds(cameras: list[Camera], exec_timeout: float = 30.0):
        if not cameras:
            return

        # background folders
        next_scan_folder, _ = cameras[0].next_scan_folder()
        bkg_target_folders: list[SysPath] = [os.path.join(next_scan_folder, f'{cam.get_name()}_Background')
                                             for cam in cameras]

        # save images
        GeecsDevice.run_no_scan(monitoring_device=cameras[0], comment=', '.join([cam.get_name() for cam in cameras])
                                + ': background collection', timeout=exec_timeout)

        # image averages
        bkg_source_paths: list[SysPath] = [os.path.join(next_scan_folder, cam.get_name()) for cam in cameras]
        for it in range(len(cameras)):
            Camera.calculate_average_image(cameras[it], bkg_source_paths[it], bkg_target_folders[it])

    @staticmethod
    def calculate_average_image(camera: Camera, bkg_imgs_folder: SysPath,
                                avg_bkg_folder: SysPath, file_name: str = 'avg_bkg.png', n_images: int = 0):
        images = glob.glob(os.path.join(bkg_imgs_folder, '*.png'))
        if n_images > 0:
            images = images[-n_images:]
        if images:
            try:
                avg_image = cv2.imread(images[0], cv2.IMREAD_GRAYSCALE)
                if len(images) > 1:
                    for it in range(len(images) - 1):
                        image_data = cv2.imread(images[it + 1], cv2.IMREAD_GRAYSCALE)
                        alpha = 1.0 / (it + 2)
                        beta = 1.0 - alpha
                        avg_image = cv2.addWeighted(image_data, alpha, avg_image, beta, 0.0)

                if not os.path.isdir(avg_bkg_folder):
                    os.makedirs(avg_bkg_folder)
                bkg_filepath: SysPath = os.path.join(avg_bkg_folder, file_name)
                cv2.imwrite(bkg_filepath, avg_image)
                camera.set(camera.var_bkg_path, value=bkg_filepath, exec_timeout=10., sync=True)

            except Exception:
                pass
