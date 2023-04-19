import queue
import re
import inspect
import time
import os
import numpy as np
from queue import Queue
import numpy.typing as npt
from threading import Thread, Condition, Event, Lock
from typing import Optional, Any, Union
from datetime import datetime as dtime
from geecs_api.api_defs import VarDict, ExpDict, VarAlias, AsyncResult, ThreadInfo
import geecs_api.interface.message_handling as mh
from geecs_api.interface import GeecsDatabase, UdpHandler, TcpSubscriber, ErrorAPI, api_error


class GeecsDevice:
    # Static
    threads_lock = Lock()
    all_threads: list[ThreadInfo] = []
    appdata_path = os.path.join(os.getenv('LOCALAPPDATA'), 'GEECS')
    scan_file_path = os.path.join(appdata_path, 'geecs_scan.txt')

    def __init__(self, name: str, exp_info: Optional[dict[str, Any]], virtual=False):

        # Static properties
        self.__dev_name: str = name.strip()
        self.__dev_virtual = virtual or not self.__dev_name
        self.__class_name = re.search(r'\w+\'>$', str(self.__class__))[0][:-2]

        # Communications
        self.mc_port: int = exp_info['MC_port']
        self.dev_tcp: Optional[TcpSubscriber] = None
        self.dev_udp: Optional[UdpHandler]
        if not self.__dev_virtual:
            self.dev_udp = UdpHandler(owner=self)
        else:
            self.dev_udp = None

        self.dev_ip: str = ''
        self.dev_port: int = 0

        # Variables
        self.dev_vars = {}
        self.var_names_by_index: dict[int, tuple[str, VarAlias]] = {}
        self.var_aliases_by_name: dict[str, tuple[VarAlias, int]] = {}

        self.setpoints: dict[VarAlias, Any] = {}
        self.state: dict[VarAlias, Any] = {}
        self.generic_vars = ['device status', 'device error', 'device preset']

        # Message handling
        self.queue_cmds = Queue()
        self.own_threads: list[(Thread, Event)] = []

        self.notify_on_udp = False
        self.queue_udp_msgs = Queue()
        self.notifier_udp_msgs = Condition()

        self.notify_on_tcp = False
        self.queue_tcp_msgs = Queue()
        self.notifier_tcp_msgs = Condition()

        if not self.__dev_virtual:
            self.dev_ip, self.dev_port = GeecsDatabase.find_device(self.__dev_name)
            if self.is_valid():
                # print(f'Device "{self.dev_name}" found: {self.dev_ip}, {self.dev_port}')
                try:
                    self.dev_tcp = TcpSubscriber(owner=self)
                    self.connect_var_listener()
                except Exception:
                    api_error.error('Failed creating TCP subscriber', 'GeecsDevice class, method "__init__"')
            else:
                api_error.warning(f'Device "{self.__dev_name}" not found', 'GeecsDevice class, method "__init__"')

            self.list_variables(exp_info['devices'])

        # Data
        self.data_root_path = exp_info['data_path']

        if not os.path.exists(GeecsDevice.appdata_path):
            os.makedirs(GeecsDevice.appdata_path)

    def cleanup(self):
        mh.flush_queue(self.queue_udp_msgs)
        mh.flush_queue(self.queue_tcp_msgs)

        self.stop_waiting_for_all_cmds()

        if self.dev_udp:
            self.dev_udp.cleanup()

        if self.dev_tcp:
            self.dev_tcp.cleanup()

    def is_valid(self):
        return not self.__dev_virtual and self.dev_ip and self.dev_port > 0

    # Accessors
    # -----------------------------------------------------------------------------------------------------------
    def get_name(self):
        return self.__dev_name

    def get_class(self):
        return self.__class_name

    # Registrations
    # -----------------------------------------------------------------------------------------------------------
    def connect_var_listener(self):
        if self.is_valid() and not self.is_var_listener_connected():
            self.dev_tcp.connect((self.dev_ip, self.dev_port))

            if not self.dev_tcp.is_connected():
                api_error.warning('Failed to connect TCP subscriber', f'GeecsDevice "{self.__dev_name}"')

    def is_var_listener_connected(self):
        return self.dev_tcp and self.dev_tcp.is_connected()

    def register_var_listener_handler(self):
        return self.dev_tcp.register_handler()  # happens only if is_valid()

    def unregister_var_listener_handler(self):
        return self.dev_tcp.unregister_handler()

    def register_cmd_executed_handler(self):
        return self.dev_udp.register_handler()  # happens only if is_valid()

    def unregister_cmd_executed_handler(self):
        return self.dev_udp.unregister_handler()

    def subscribe_var_values(self, variables: Optional[list[str]] = None) -> bool:
        subscribed = False

        if self.is_valid() and variables is None:
            variables = [var[0] for var in self.var_names_by_index.values()]

        variables = self.generic_vars + variables

        if self.is_valid() and variables:
            try:
                subscribed = self.dev_tcp.subscribe(','.join(variables))
            except Exception as ex:
                api_error.error(str(ex), 'Class GeecsDevice, method "subscribe_var_values"')

        return subscribed

    def unsubscribe_var_values(self):
        if self.is_var_listener_connected():
            self.dev_tcp.unsubscribe()

    # Variables
    # -----------------------------------------------------------------------------------------------------------
    def list_variables(self, exp_devs: Optional[ExpDict] = None) -> tuple[VarDict, ExpDict]:
        try:
            if exp_devs is None:
                exp_info = GeecsDatabase.collect_exp_info()
                exp_devs = exp_info['devices']

            self.dev_vars = exp_devs[self.__dev_name]

        except Exception:
            self.dev_vars = {}

        return self.dev_vars, exp_devs

    def find_var_by_alias(self, alias: VarAlias = VarAlias('')) -> str:
        if not self.dev_vars:
            self.list_variables()

        if not self.dev_vars:
            return ''

        var_name = ''
        for attributes in self.dev_vars.values():
            if attributes['alias'] == alias:
                var_name = attributes['variablename']
                break

        if not var_name and alias in self.dev_vars:
            var_name = str(alias)

        return var_name

    def build_var_dicts(self, aliases: tuple[VarAlias]):
        self.var_names_by_index: dict[int, tuple[str, VarAlias]] = \
            {index: (self.find_var_by_alias(aliases[index]), aliases[index]) for index in range(len(aliases))}

        self.var_aliases_by_name: dict[str, tuple[VarAlias, int]] = \
            {self.find_var_by_alias(aliases[index]): (aliases[index], index) for index in range(len(aliases))}

    def state_value(self, var_name: str) -> Any:
        var_alias: VarAlias

        if var_name in self.generic_vars:
            var_alias = VarAlias(var_name)
        else:
            var_alias = self.var_aliases_by_name[var_name][0]

        if var_alias in self.state:
            return self.state[var_alias]
        else:
            return None

    # Operations
    # -----------------------------------------------------------------------------------------------------------
    def set(self, variable: str, value, exec_timeout: Optional[float] = 120.0, attempts_max: int = 5, sync=True)\
            -> AsyncResult:
        return self._execute(variable, value, exec_timeout, attempts_max, sync)

    def get(self, variable: str, exec_timeout: Optional[float] = 5.0, attempts_max: int = 5, sync=True)\
            -> AsyncResult:
        return self._execute(variable, None, exec_timeout, attempts_max, sync)

    def _start_scan(self, timeout: float = 300.) -> tuple[bool, bool]:
        cmd = f'FileScan>>{GeecsDevice.scan_file_path}'
        accepted = self.dev_udp.send_scan_cmd(cmd)

        time.sleep(20.)  # to enter info (won't be needed in the future, hopefully) and devices to enter scan mode
        t0 = time.monotonic()
        while True:
            timed_out = (time.monotonic() - t0 > timeout)
            if (self.state[VarAlias('device status')] == 'scan') or timed_out:
                break
            time.sleep(1.)

        return accepted, timed_out

    def get_status(self, exec_timeout: float = 2.0, sync=True) -> Union[Optional[float], AsyncResult]:
        ret = self.get('device status', exec_timeout=exec_timeout, sync=sync)
        if sync:
            return self.state_value('device status')
        else:
            return ret

    def interpret_value(self, var_alias: VarAlias, val_string: str) -> Any:
        return float(val_string)

    def interpret_generic_variables(self, var: str, val: str):
        # ['device status', 'device error', 'device preset']
        self.state[VarAlias(var)] = val

    def _execute(self, variable: str, value, exec_timeout: Optional[float] = 10.0,
                 attempts_max: int = 5, sync=True) -> AsyncResult:
        if api_error.is_error:
            return False, '', (None, None)

        scan = (variable == 'scan')

        if scan:
            cmd_str = f'StartScan>>{GeecsDevice.scan_file_path}'
            cmd_label = 'scan'
        else:
            if isinstance(value, float):
                cmd_str = f'set{variable}>>{value:.6f}'
                cmd_label = f'set({variable}, {value:.6f})'
            elif isinstance(value, str):
                cmd_str = f'set{variable}>>{value}'
                cmd_label = f'set({variable}, {value})'
            elif isinstance(value, bool):
                cmd_str = f'set{variable}>>{int(value)}'
                cmd_label = f'set({variable}, {value})'
            else:
                cmd_str = f'get{variable}>>'
                cmd_label = f'get({variable})'

        if not self.is_valid():
            api_error.warning(f'Failed to execute "{cmd_label}"',
                              f'GeecsDevice "{self.__dev_name}" not connected')
            return False, '', (None, None)

        stamp = re.sub(r'[\s.:]', '-', dtime.now().__str__())
        cmd_label += f' @ {stamp}'

        queued: bool = False
        async_thread: ThreadInfo = (None, None)

        self._cleanup_threads()

        if sync or scan:
            self.wait_for_all_cmds(timeout=120.)

            with GeecsDevice.threads_lock:
                self.process_command(cmd_str, cmd_label, thread_info=(None, None), attempts_max=attempts_max)
                if not scan:
                    self.dev_udp.cmd_checker.wait_for_exe(cmd_tag=cmd_label, timeout=exec_timeout, sync=sync)

        elif exec_timeout > 0:
            with GeecsDevice.threads_lock:
                # create listening thread (only)
                async_thread: ThreadInfo = \
                    self.dev_udp.cmd_checker.wait_for_exe(cmd_tag=cmd_label, timeout=exec_timeout, sync=sync)

                # if nothing running and no commands in queue
                if (not self.own_threads) and self.queue_cmds.empty():
                    self.process_command(cmd_str, cmd_label, thread_info=async_thread, attempts_max=attempts_max)
                else:
                    self.queue_cmds.put((cmd_str, cmd_label, async_thread, attempts_max))
                    queued = True

        return queued, cmd_label, async_thread

    def dequeue_command(self):
        self._cleanup_threads()

        with GeecsDevice.threads_lock:
            # if nothing running and commands in queue
            if (not self.own_threads) and (not self.queue_cmds.empty()):
                try:
                    cmd_str, cmd_label, async_thread, attempts_max = self.queue_cmds.get_nowait()
                    self.process_command(cmd_str, cmd_label, thread_info=async_thread, attempts_max=attempts_max)
                except queue.Empty:
                    pass

    def process_command(self, cmd_str: str, cmd_label: str,
                        thread_info: ThreadInfo = (None, None), attempts_max: int = 5):
        accepted = False
        try:
            for _ in range(attempts_max):
                sent = self.dev_udp.send_cmd(ipv4=(self.dev_ip, self.dev_port), msg=cmd_str)
                if sent:
                    accepted = self.dev_udp.ack_cmd(timeout=5.0)
                else:
                    time.sleep(0.1)
                    continue

                if accepted or api_error.is_error:
                    break
                else:
                    time.sleep(0.1)

        except Exception as ex:
            api_error.error(str(ex), f'GeecsDevice "{self.__dev_name}", method "{cmd_label}"')

        if accepted and (thread_info[0] is not None):
            thread_info[0].start()
            self.own_threads.append(thread_info)
            GeecsDevice.all_threads.append(thread_info)

    def handle_response(self, net_msg: mh.NetworkMessage) -> tuple[str, str, str, bool]:
        try:
            dev_name, cmd_received, dev_val, err_status = GeecsDevice._response_parser(net_msg.msg)

            # Queue & notify
            if self.notify_on_udp:
                self.queue_udp_msgs.put((dev_name, cmd_received, dev_val, err_status))
                self.notifier_udp_msgs.notify_all()

            # Error handling
            if net_msg.err.is_error or net_msg.err.is_warning:
                print(net_msg.err)

            if dev_name != self.__dev_name:
                warn = ErrorAPI('Mismatch in device name', f'Class {self.__class_name}, method "handle_response"')
                print(warn)

            # Update dictionaries
            if dev_name == self.get_name() and not err_status and cmd_received[:3] == 'get':
                var_alias = VarAlias('')

                if cmd_received[3:] in self.generic_vars:
                    self.interpret_generic_variables(cmd_received[3:], dev_val)
                    var_alias = VarAlias(cmd_received[3:])

                if cmd_received[3:] in self.var_aliases_by_name:
                    var_alias = self.var_aliases_by_name[cmd_received[3:]][0]
                    dev_val = self.interpret_value(var_alias, dev_val)
                    self.state[var_alias] = dev_val

                if var_alias:
                    dev_val = f'"{dev_val}"' if isinstance(dev_val, str) else dev_val
                    print(f'{self.__class_name} [{self.__dev_name}]: {var_alias} = {dev_val}')

            if dev_name == self.get_name() and not err_status and cmd_received[:3] == 'set':
                var_alias = VarAlias('')

                if cmd_received[3:] in self.var_aliases_by_name:
                    var_alias = self.var_aliases_by_name[cmd_received[3:]][0]
                    dev_val = self.interpret_value(var_alias, dev_val)
                    self.setpoints[var_alias] = dev_val

                if var_alias:
                    dev_val = f'"{dev_val}"' if isinstance(dev_val, str) else dev_val
                    print(f'{self.__class_name} [{self.__dev_name}]: {var_alias} set to {dev_val}')

            return dev_name, cmd_received, dev_val, err_status

        except Exception as ex:
            err = ErrorAPI(str(ex), f'Class {self.__class_name}, method "{inspect.stack()[0][3]}"')
            print(err)
            return '', '', '', True

    def handle_subscription(self, net_msg: mh.NetworkMessage) -> tuple[str, int, dict[str, str]]:
        try:
            dev_name, shot_nb, dict_vals = GeecsDevice._subscription_parser(net_msg.msg)

            # Queue & notify
            if self.notify_on_tcp:
                self.queue_tcp_msgs.put((dev_name, shot_nb, dict_vals))
                self.notifier_tcp_msgs.notify_all()

            # Error handling
            if net_msg.err.is_error or net_msg.err.is_warning:
                print(net_msg.err)

            # Update dictionaries
            if dev_name == self.get_name() and dict_vals:
                for var, val in dict_vals.items():
                    if var in self.generic_vars:
                        self.interpret_generic_variables(var, val)

                    if var in self.var_aliases_by_name:
                        var_alias: VarAlias = self.var_aliases_by_name[var][0]
                        self.state[var_alias] = self.interpret_value(var_alias, val)

            return dev_name, shot_nb, dict_vals

        except Exception as ex:
            err = ErrorAPI(str(ex), f'Class {self.__class_name}, method "{inspect.stack()[0][3]}"')
            print(err)
            return '', 0, {}

    @staticmethod
    def _subscription_parser(msg: str = '') -> tuple[str, int, dict[str, str]]:
        """ General parser to be called when messages are received. """

        # msg = 'U_S2V>>0>>Current nval, -0.000080 nvar, Voltage nval,0.002420 nvar,'
        pattern = re.compile(r'[^,]+nval,[^,]+nvar')
        blocks = msg.split('>>')
        dev_name = blocks[0]
        shot_nb = int(blocks[1])
        vars_vals = pattern.findall(blocks[-1])

        dict_vals = {vars_vals[i].split(',')[0][:-5].strip(): vars_vals[i].split(',')[1][:-5]
                     for i in range(len(vars_vals))}

        return dev_name, shot_nb, dict_vals

    @staticmethod
    def _response_parser(msg: str = '') -> tuple[str, str, str, bool]:
        """ General parser to be called when messages are received. """

        # Examples:
        # 'U_ESP_JetXYZ>>getJet_X (mm)>>>>error,Error occurred during access CVT -  "jet_x (mm)" variable not found'
        # 'U_ESP_JetXYZ>>getPosition.Axis 1>>7.600390>>no error,'

        dev_name, cmd_received, dev_val, err_msg = msg.split('>>')
        err_status, err_msg = err_msg.split(',')
        err_status = (err_status == 'error')
        if err_status:
            api_error.error(err_msg, f'Failed to execute command "{cmd_received}", error originated in control system')

        return dev_name, cmd_received, dev_val, err_status

    def coerce_float(self, var_alias: VarAlias, method: str, value: float,
                     span: tuple[Optional[float], Optional[float]]) -> float:
        try:
            if span[0] and value < span[0]:
                api_error.warning(f'{var_alias} value coerced from {value} to {span[0]}',
                                  f'Class {self.__class_name}, method "{method}"')
                value = span[0]
            if span[1] and value > span[1]:
                api_error.warning(f'{var_alias} value coerced from {value} to {span[1]}',
                                  f'Class {self.__class_name}, method "{method}"')
                value = span[1]
        except Exception:
            api_error.error('Failed to coerce value')

        return value

    def _scan_values(self, var_alias: VarAlias, start_value: float, end_value: float, step_size: float,
                     spans:  dict[VarAlias, tuple[float, float]]) -> npt.ArrayLike:
        start_value = self.coerce_float(var_alias, inspect.stack()[0][3], start_value, spans[var_alias])
        end_value = self.coerce_float(var_alias, inspect.stack()[0][3], end_value, spans[var_alias])
        if end_value < start_value:
            step_size = -abs(step_size)
        else:
            step_size = abs(step_size)
        return np.arange(start_value, end_value + step_size, step_size)

    @staticmethod
    def _write_scan_file(devices: Union[list[str], str], variables: Union[list[str], str],
                         values_by_row: npt.ArrayLike, shots_per_step: int = 10):
        scan_number = 1

        with open(GeecsDevice.scan_file_path, 'w+') as f:
            f.write(f'[Scan{scan_number}]\n')
            if isinstance(devices, list):
                f.write('Device = "' + ','.join(devices) + '"\n')
            else:
                f.write(f'Device = "{devices}"\n')

            if isinstance(variables, list):
                f.write('Variable = "' + ','.join(variables) + '"\n')
            else:
                f.write(f'Variable = "{variables}"\n')
            f.write('Values:#shots = "')

            if values_by_row.ndim > 1:
                for col in range(values_by_row.shape[1]):
                    f.write(f'({str(list(values_by_row[:, col]))[1:-1]}):{shots_per_step}|')
            else:
                for col in range(values_by_row.size):
                    f.write(f'({values_by_row[col]}):{shots_per_step}|')
            f.write('"')

    # Synchronization
    # -----------------------------------------------------------------------------------------------------------
    @staticmethod
    def cleanup_all_threads():
        with GeecsDevice.threads_lock:
            for it in range(len(GeecsDevice.all_threads)):
                if not GeecsDevice.all_threads[-1 - it][0].is_alive():
                    GeecsDevice.all_threads.pop(-1 - it)

    def _cleanup_threads(self):
        with GeecsDevice.threads_lock:
            for it in range(len(self.own_threads)):
                try:
                    if not self.own_threads[-1 - it][0].is_alive():
                        self.own_threads.pop(-1 - it)
                except Exception:
                    continue

            for it in range(len(GeecsDevice.all_threads)):
                try:
                    if not GeecsDevice.all_threads[-1 - it][0].is_alive():
                        GeecsDevice.all_threads.pop(-1 - it)
                except Exception:
                    continue

    @staticmethod
    def wait_for_all_devices(timeout: Optional[float] = None):
        GeecsDevice.cleanup_all_threads()
        synced = True

        with GeecsDevice.threads_lock:
            for thread in GeecsDevice.all_threads:
                thread[0].join(timeout)
                synced &= thread[0].is_alive()

        return synced

    @staticmethod
    def stop_waiting_for_all_devices():
        GeecsDevice.cleanup_all_threads()

        with GeecsDevice.threads_lock:
            for thread in GeecsDevice.all_threads:
                thread[1].set()

    def wait_for_all_cmds(self, timeout: Optional[float] = None) -> bool:
        self._cleanup_threads()
        any_alive = False

        with GeecsDevice.threads_lock:
            for thread in self.own_threads:
                thread[0].join(timeout)
                any_alive |= thread[0].is_alive()

        return not any_alive

    def stop_waiting_for_all_cmds(self):
        self._cleanup_threads()

        with GeecsDevice.threads_lock:
            for thread in self.own_threads:
                thread[1].set()

    def wait_for_cmd(self, thread: Thread, timeout: Optional[float] = None):
        with GeecsDevice.threads_lock:
            if self.is_valid() and thread.is_alive():
                thread.join(timeout)

            alive = thread.is_alive()

        return not alive

    def stop_waiting_for_cmd(self, thread: Thread, stop: Event):
        if self.is_valid() and thread.is_alive():
            stop.set()
