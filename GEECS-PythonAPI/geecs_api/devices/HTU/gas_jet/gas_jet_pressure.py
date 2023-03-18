from __future__ import annotations
import time
from queue import Queue
from typing import Optional, Any
from threading import Thread, Condition, Event
import geecs_api.interface.message_handling as mh
from geecs_api.devices.geecs_device import GeecsDevice
from geecs_api.interface import GeecsDatabase, ErrorAPI, api_error


class GasJetPressure(GeecsDevice):
    # Singleton
    __instance = None

    def __new__(cls, *args, **kwargs):
        if cls.__instance is None:
            cls.__instance = super(GasJetPressure, cls).__new__(cls)
            cls.__instance.__initialized = False
        return cls.__instance

    def __init__(self, exp_vars: dict[str, dict[str, dict[str, Any]]]):
        if self.__initialized:
            return
        self.__initialized = True

        super().__init__('U_HP_Daq')

        self.list_variables(exp_vars)
        aliases = ['',
                   '']
        self.get_var_dicts(aliases)

        self.register_cmd_executed_handler()
        self.register_var_listener_handler()

    def handle_response(self, net_msg: mh.NetworkMessage,
                        notifier: Optional[Condition] = None,
                        queue_msgs: Optional[Queue] = None) -> tuple[str, str, str, str]:
        try:
            dev_name, cmd_received, dev_val, err_status = super().handle_response(net_msg, notifier, queue_msgs)

            if dev_name == self.dev_name and dev_val and cmd_received[:3] == 'get':
                var_alias = self.var_aliases[cmd_received[3:]][0]
                self.gets[var_alias] = float(dev_val)
                print(f'{var_alias} = {float(dev_val)}')

            if dev_name == self.dev_name and dev_val and cmd_received[:3] == 'set':
                var_alias = self.var_aliases[cmd_received[3:]][0]
                self.sets[var_alias] = float(dev_val)
                print(f'{var_alias} set to {float(dev_val)}')

            return dev_name, cmd_received, dev_val, err_status

        except Exception as ex:
            err = ErrorAPI(str(ex), f'Class {self.__class__}, method "handle_response"')
            print(err)

    def handle_subscription(self, net_msg: mh.NetworkMessage,
                            notifier: Optional[Condition] = None,
                            queue_msgs: Optional[Queue] = None) -> tuple[str, int, dict[str, float]]:
        try:
            dev_name, shot_nb, dict_vals = super().handle_subscription(net_msg, notifier, queue_msgs)

            if dev_name == self.dev_name and dict_vals:
                for var, val in dict_vals.items():
                    if var in self.var_aliases:
                        var_alias = self.var_aliases[var][0]
                        self.gets[var_alias] = float(val)

            return dev_name, shot_nb, dict_vals

        except Exception as ex:
            err = ErrorAPI(str(ex), 'Class GeecsDevice, method "subscription_handler"')
            print(err)
            return '', 0, {}

    def get_axis_var_name(self, axis: int):
        if axis < 0 or axis > 2:
            return ''
        else:
            return self.var_names.get(axis)[0]

    def get_axis_var_alias(self, axis: int):
        if axis < 0 or axis > 2:
            return ''
        else:
            return self.var_names.get(axis)[1]

    def get_position(self, axis: Optional[str, int], exec_timeout: float = 2.0, sync=True) \
            -> tuple[bool, str, tuple[Optional[Thread], Optional[Event]]]:
        if isinstance(axis, str):
            if len(axis) == 1:
                axis = ord(axis.upper()) - ord('X')
            else:
                axis = -1

        if axis < 0 or axis > 2:
            return False, '', (None, None)

        return self.get(self.get_axis_var_name(axis), exec_timeout=exec_timeout, sync=sync)

    def set_position(self, axis: Optional[str, int], value: float, exec_timeout: float = 2.0, sync=True) \
            -> tuple[bool, str, tuple[Optional[Thread], Optional[Event]]]:
        if isinstance(axis, str):
            if len(axis) == 1:
                axis = ord(axis.upper()) - ord('X')
            else:
                axis = -1

        if axis < 0 or axis > 2:
            return False, '', (None, None)

        return self.set(self.get_axis_var_name(axis), value=value, exec_timeout=exec_timeout, sync=sync)


if __name__ == '__main__':
    api_error.clear()

    # list experiment devices and variables
    exp_devs = GeecsDatabase.find_experiment_variables('Undulator')

    # create gas jet object
    jet_pressure = GasJetPressure(exp_devs)
    print(f'Variables subscription: {jet_pressure.subscribe_var_values()}')

    # get X-position
    time.sleep(1.0)
    # jet.get_position('X', sync=True)

    # retrieve currently known positions
    try:
        print(f'Pressure state:\n\t{jet_pressure.gets}')
        print(f'Pressure config:\n\t{jet_pressure.sets}')
    except Exception as e:
        api_error.error(str(e), 'Demo code for gas jet')
        pass

    jet_pressure.cleanup()
    print(api_error)
