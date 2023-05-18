from __future__ import annotations
import os
import time
import shutil
import cv2
import numpy as np
from typing import Any
from datetime import datetime as dtime
from geecs_api.tools.images.batches import average_images
from geecs_api.api_defs import VarAlias, SysPath
from geecs_api.devices.geecs_device import GeecsDevice


class Camera(GeecsDevice):
    # ROIs with [left, right, top, bottom] (x_lim = [:1], y_lim = [-2:])
    ROIs = {'UC_ALineEbeam1': [204, 777, 319, 701],
            'UC_ALineEBeam2': [261, 843, 499, 656],
            'UC_ALineEBeam3': [274, 528, 256, 546],
            'UC_VisaEBeam1': [105, 708, 252, 377],
            'UC_VisaEBeam2': [164, 434, 100, 400],
            'UC_VisaEBeam3': [185, 477, 137, 469],
            'UC_VisaEBeam4': [263, 541, 192, 508],
            'UC_VisaEBeam5': [167, 450, 106, 427],
            'UC_VisaEBeam6': [167, 500, 147, 406],
            'UC_VisaEBeam7': [125, 462, 130, 490],
            'UC_VisaEBeam8': [206, 401, 111, 466],
            'UC_VisaEBeam9': [341, 1141, 628, 670],
            'UC_UndulatorRad2': [600, 2360, 1420, 1170]}

    def __init__(self, device_name: str):
        super().__init__(device_name)

        # self.gui_path: SysPath = GeecsDevice.exp_info['GUIs'][device_name]

        self.bkg_folder: SysPath = os.path.join(self.data_root_path, 'backgrounds', f'{device_name}')
        if not os.path.isdir(self.bkg_folder):
            os.makedirs(self.bkg_folder)

        self.__variables = {VarAlias('BackgroundPath'): (None, None),
                            VarAlias('localsavingpath'): (None, None),
                            VarAlias('exposure'): (None, None),
                            VarAlias('triggerdelay'): (None, None)}
        self.build_var_dicts(tuple(self.__variables.keys()))
        self.var_bkg_path: str = self.var_names_by_index.get(0)[0]
        self.var_save_path: str = self.var_names_by_index.get(1)[0]
        self.var_exposure: str = self.var_names_by_index.get(2)[0]

        if self.get_name() in Camera.ROIs:
            self.roi = np.array(Camera.ROIs[self.get_name()])
        else:
            self.roi = None

        self.label = Camera.label_from_name(self.get_name())
        if self.label in ['A1', 'A2']:
            self.rot_90: int = 90
        else:
            self.rot_90 = 0

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

    def save_local_background(self, n_images: int = 15, set_as_background: bool = True):
        # background images folder
        source_path: SysPath = os.path.join(self.bkg_folder, 'tmp_images')
        if os.path.isdir(source_path):
            shutil.rmtree(source_path, ignore_errors=True)
        while True:
            if not os.path.isdir(source_path):
                break
        os.makedirs(source_path)

        self.set(self.var_save_path, value=source_path, exec_timeout=10., sync=True)

        # file name
        stamp = dtime.now()
        file_name = stamp.strftime('%y%m%d_%H%M') + f'_{self.get_name()}_x{n_images}_Local.png'
        file_path: SysPath = os.path.join(self.bkg_folder, file_name)

        # save images
        if n_images > 0:
            self.set('save', 'on', exec_timeout=5.)
            time.sleep(n_images + 2)
            self.set('save', 'off', exec_timeout=5.)

            # wait to write files to disk
            t0 = time.monotonic()
            while True:
                if time.monotonic() - t0 > 10. or len(next(os.walk(source_path))[2]) >= n_images:
                    break

        # average image
        avg_image, _ = average_images(source_path, n_images)

        if avg_image:
            cv2.imwrite(file_path, avg_image)
            if set_as_background:
                time.sleep(1.)  # buffer to write file to disk
                self.set(self.var_bkg_path, value=file_path, exec_timeout=10., sync=True)

    def save_background(self, exec_timeout: float = 30.0, set_as_background: bool = True):
        next_scan_folder, _ = self.next_scan_folder()
        scan_name = os.path.basename(next_scan_folder)

        # background file name
        stamp = dtime.now()
        file_name: str = stamp.strftime('%y%m%d_%H%M') + f'_{self.get_name()}_{scan_name}.png'
        file_path: SysPath = os.path.join(self.bkg_folder, file_name)

        # save images
        GeecsDevice.run_no_scan(monitoring_device=self,
                                comment=f'{self.get_name()}: background collection',
                                timeout=exec_timeout)

        # average image
        avg_image, _ = average_images(images_folder=os.path.join(next_scan_folder, self.get_name()))

        if avg_image:
            cv2.imwrite(file_path, avg_image)
            if set_as_background:
                time.sleep(1.)  # buffer to write file to disk
                self.set(self.var_bkg_path, value=file_path, exec_timeout=10., sync=True)

    @staticmethod
    def label_from_name(name: str):
        if name[3] == 'A':
            return f'A{name[-1]}'
        elif name[3] == 'V':
            return f'U{name[-1]}'
        elif name[3] == 'U':
            return name[-4:]
        elif name[3] == 'D':
            return 'DC'
        elif name[3] == 'P':
            return 'P1'
        else:
            return name

    @staticmethod
    def name_from_label(label: str):
        labels = [Camera.label_from_name(name) for name in Camera.ROIs.keys()]
        if label in labels:
            return list(Camera.ROIs.keys())[labels.index(label)]
        else:
            return label

    @staticmethod
    def save_multiple_backgrounds(cameras: list[Camera], exec_timeout: float = 30.0, set_as_background: bool = True):
        if not cameras:
            return
        next_scan_folder, _ = cameras[0].next_scan_folder()

        # background file name
        stamp = dtime.now()
        file_stamp = stamp.strftime('%y%m%d_%H%M')
        scan_name = os.path.basename(next_scan_folder)

        # save images
        GeecsDevice.run_no_scan(monitoring_device=cameras[0],
                                comment=', '.join([cam.get_name() for cam in cameras]) + ': background collection',
                                timeout=exec_timeout)

        # average images
        for camera in cameras:
            try:
                file_name = f'{file_stamp}_{camera.get_name()}_{scan_name}.png'
                file_path: SysPath = os.path.join(camera.bkg_folder, file_name)

                avg_image, _ = average_images(images_folder=os.path.join(next_scan_folder, camera.get_name()))

                if avg_image:
                    cv2.imwrite(file_path, avg_image)
                    if set_as_background:
                        time.sleep(1.)  # buffer to write file to disk
                        camera.set(camera.var_bkg_path, value=file_path, exec_timeout=10., sync=True)

            except Exception:
                continue
