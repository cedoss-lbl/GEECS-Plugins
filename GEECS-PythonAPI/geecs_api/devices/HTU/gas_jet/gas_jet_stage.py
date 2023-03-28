from __future__ import annotations
import time
import inspect
from typing import Optional, Any
from threading import Thread, Event
from geecs_api.api_defs import *
from geecs_api.devices.geecs_device import GeecsDevice
from geecs_api.interface import GeecsDatabase, api_error


class GasJetStage(GeecsDevice):
    # Singleton
    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, 'instance'):
            cls.instance = super(GasJetStage, cls).__new__(cls)
            cls.instance.__initialized = False
        return cls.instance

    def __init__(self, exp_vars: dict[str, dict[str, dict[str, Any]]]):
        if self.__initialized:
            return
        self.__initialized = True

        super().__init__('U_ESP_JetXYZ', exp_vars)

        self.__variables = {VarAlias('Jet_X (mm)'): (2., 10.),  # [min, max]
                            VarAlias('Jet_Y (mm)'): (-8., -1.),
                            VarAlias('Jet_Z (mm)'): (0., 25.)}
        self.build_var_dicts(tuple(self.__variables.keys()))

        self.register_cmd_executed_handler()
        self.register_var_listener_handler()

    def get_axis_var_name(self, axis: int) -> str:
        if axis < 0 or axis > 2:
            return ''
        else:
            return self.var_names_by_index.get(axis)[0]

    def get_axis_var_alias(self, axis: int) -> VarAlias:
        if axis < 0 or axis > 2:
            return VarAlias('')
        else:
            return self.var_names_by_index.get(axis)[1]

    def state_x(self) -> Optional[float]:
        return self.state_value(self.get_axis_var_name(0))

    def state_y(self) -> Optional[float]:
        return self.state_value(self.get_axis_var_name(1))

    def state_z(self) -> Optional[float]:
        return self.state_value(self.get_axis_var_name(2))

    def get_position(self, axis: Optional[str, int], exec_timeout: float = 2.0, sync=True) \
            -> Union[Optional[float], AsyncResult]:
        if isinstance(axis, str):
            if len(axis) == 1:
                axis = ord(axis.upper()) - ord('X')
            else:
                axis = -1

        if axis < 0 or axis > 2:
            if sync:
                return None
            else:
                return False, '', (None, None)

        ret = self.get(self.get_axis_var_name(axis), exec_timeout=exec_timeout, sync=sync)
        if sync:
            return self.state_value(self.get_axis_var_name(axis))
        else:
            return ret

    def set_position(self, axis: Optional[str, int], value: float, exec_timeout: float = 30.0, sync=True) \
            -> Union[Optional[float], AsyncResult]:
        if isinstance(axis, str):
            if len(axis) == 1:
                axis = ord(axis.upper()) - ord('X')
            else:
                axis = -1

        if axis < 0 or axis > 2:
            if sync:
                return None
            else:
                return False, '', (None, None)

        var_alias = self.get_axis_var_alias(axis)
        value = self.coerce_float(var_alias, inspect.stack()[0][3], value, self.__variables[var_alias])

        ret = self.set(self.get_axis_var_name(axis), value=value, exec_timeout=exec_timeout, sync=sync)
        if sync:
            return self.state_value(self.get_axis_var_name(axis))
        else:
            return ret


if __name__ == '__main__':
    api_error.clear()

    # list experiment devices and variables
    exp_devs = GeecsDatabase.find_experiment_variables('Undulator')

    # create gas jet object
    jet = GasJetStage(exp_devs)
    other_jet = GasJetStage(exp_devs)
    print(f'Only one jet: {jet is other_jet}')
    print(f'Variables subscription: {jet.subscribe_var_values()}')

    # set position
    time.sleep(1.0)
    print(f'Jet state: {jet.state}')
    # jet.get_position('Y', sync=True)
    jet.set_position('Y', jet.state[jet.get_axis_var_alias(1)])

    print(f'Jet state: {jet.state}')
    jet.cleanup()
    print(api_error)