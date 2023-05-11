from __future__ import annotations
from geecs_api.api_defs import VarAlias
from geecs_api.devices.geecs_device import GeecsDevice
from geecs_api.interface import GeecsDatabase, api_error
from geecs_api.devices.HTU.multi_channels import PlungersPLC
from geecs_api.devices.HTU.diagnostics.cameras import Camera
from geecs_api.devices.HTU.diagnostics.ebeam_phosphor.phosphors import Phosphor


class EBeamPhosphor(GeecsDevice):
    """ e-beam diagnostic made of a camera-phosphors plunger pair. """
    def __init__(self, camera_name: str, plunger_controller: GeecsDevice, plunger_name: str,
                 tcp_subscription: bool = True):
        super().__init__('e_beam_phosphor', virtual=True)

        self.camera = Camera(camera_name)
        self.screen = Phosphor(f'{plunger_controller.get_name()}_{plunger_name}',
                               VarAlias(plunger_name), plunger_controller)

        if tcp_subscription:
            self.camera.subscribe_var_values()
            self.screen.subscribe_var_values()

    def close(self):
        self.camera.close()
        self.screen.close()


if __name__ == '__main__':
    api_error.clear()
    GeecsDevice.exp_info = GeecsDatabase.collect_exp_info('Undulator')

    PLC = PlungersPLC()
    e_beam_phosphor_A1 = EBeamPhosphor(camera_name='UC_ALineEbeam1',
                                       plunger_controller=PLC,
                                       plunger_name='ALine1 plunger',
                                       tcp_subscription=False)

<<<<<<< Updated upstream:GEECS-PythonAPI/geecs_api/devices/HTU/diagnostics/ebeam_phosphor/ebeam_phosphor.py
    e_beam_phosphor_A1.cleanup()
    PLC.cleanup()
=======
    e_imager_A1.close()
    PLC.close()
>>>>>>> Stashed changes:GEECS-PythonAPI/geecs_api/devices/HTU/diagnostics/e_imager/e_imager.py
