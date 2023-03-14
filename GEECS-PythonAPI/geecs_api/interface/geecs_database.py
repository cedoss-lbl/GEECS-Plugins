import os
from typing import Any
import configparser
import mysql.connector
from geecs_api.interface.geecs_errors import api_error
import tkinter as tk
from tkinter import filedialog

# Pop-ups initialization
tk_root = tk.Tk()
tk_root.withdraw()


def find_database():
    default_path = r'C:\GEECS\user data'
    default_name = 'Configurations.INI'

    db_name = db_ip = db_user = db_pwd = ''

    if not os.path.isfile(os.path.join(default_path, default_name)):
        path_cfg = filedialog.askopenfilename(filetypes=[('INI Files', '*.INI'), ('All Files', '*.*')],
                                              initialdir=default_path,
                                              initialfile=default_name,
                                              title='Choose a configuration file:')
    else:
        path_cfg = os.path.join(default_path, default_name)

    if path_cfg:
        try:
            config = configparser.ConfigParser()
            config.read(path_cfg)

            db_name = config['Database']['name']
            db_ip = config['Database']['ipaddress']
            db_user = config['Database']['user']
            db_pwd = config['Database']['password']

        except Exception:
            pass

    return db_name, db_ip, db_user, db_pwd


class GeecsDatabase:
    name, ipv4, username, password = find_database()

    @staticmethod
    def find_experiment_variables(exp_name: str = 'Undulator') -> dict[str, dict[str, dict[str, Any]]]:
        """ Dictionary of (key) devices with (values) dictionary of (key) variables and (values) attributes. """

        db = mysql.connector.connect(
            host=GeecsDatabase.ipv4,
            user=GeecsDatabase.username,
            password=GeecsDatabase.password,
            database=GeecsDatabase.name)

        db_cursor = db.cursor(dictionary=True)
        cmd_str = """
            SELECT * FROM
                -- subquery that returns devicename, variablename, and source table where the defaultvalue,
                -- min, max, etc. should come from
                (
                    SELECT devicename, variablename, MAX(precedence_sourcetable) AS precedence_sourcetable
                    FROM
                    (
                        (
                        SELECT `name` AS variablename, device AS devicename,
                        '2_variable' AS precedence_sourcetable FROM variable
                        )
                    UNION
                        (
                        SELECT devicetype_variable.name AS variablename, device.name AS devicename,
                        '1_devicetype_variable' AS precedence_sourcetable
                        FROM devicetype_variable JOIN device ON devicetype_variable.devicetype = device.devicetype
                        )
                    ) AS variable_device_from_both_tables GROUP BY devicename, variablename
                ) AS max_precedence
                -- subquery containing defaultvalue, min, max, etc from both tables. the correct version,
                -- by precedence, is selected through the join.
                LEFT JOIN
                (
                    (
                    SELECT variable.name AS variablename, variable.device AS devicename,
                    '2_variable' AS precedence_sourcetable, defaultvalue, `min`, `max`, stepsize, units,
                    choice_id, tolerance, alias, default_experiment
                    FROM variable JOIN device ON variable.device = device.name -- to pull default_experiment
                    )
                UNION
                    (
                    SELECT devicetype_variable.name AS variablename, device.name AS devicename,
                    '1_devicetype_variable' AS precedence_sourcetable, defaultvalue, `min`, `max`, stepsize, units,
                    choice_id, tolerance, alias, default_experiment
                    FROM devicetype_variable JOIN device ON devicetype_variable.devicetype = device.devicetype
                    )
                ) AS variable_device_parameters_from_both_tables 
                USING (variablename, devicename, precedence_sourcetable) 
                -- Get datatype
                LEFT JOIN (SELECT id AS choice_id, choices FROM choice) AS datatype USING (choice_id)
                -- Now filter for device, experiment
            WHERE default_experiment = %s;
        """

        db_cursor.execute(cmd_str, (exp_name,))
        rows = db_cursor.fetchall()

        exp_vars: dict[str, dict[str, dict[str, Any]]] = {}
        while rows:
            row = rows.pop()
            if row['devicename'] in exp_vars:
                exp_vars[row['devicename']][row['variablename']] = row
            else:
                exp_vars[row['devicename']] = {row['variablename']: row}

        return exp_vars

    @staticmethod
    def find_device(dev_name=''):
        db_cursor = None
        dev_ip: str = ''
        dev_port: int = 0

        try:
            db = mysql.connector.connect(
                host=GeecsDatabase.ipv4,
                user=GeecsDatabase.username,
                password=GeecsDatabase.password)

            selectors = ["ipaddress", "commport"]

            db_cursor = db.cursor()
            db_cursor.execute(f'SELECT {",".join(selectors)} FROM {GeecsDatabase.name}.device WHERE name=%s;',
                              (dev_name,))
            db_result = db_cursor.fetchone()
            dev_ip = db_result[0]
            dev_port = int(db_result[1])

        except Exception as ex:
            api_error.error(str(ex), f'GeecsDatabase class, static method "find_device({dev_name})"')

        finally:
            try:
                db_cursor.close()
            except Exception:
                pass

        return dev_ip, dev_port


if __name__ == '__main__':
    print('Name:\n\t' + GeecsDatabase.name)
    print('IP:\n\t' + GeecsDatabase.ipv4)
    print('User:\n\t' + GeecsDatabase.username)
    print('Password:\n\t' + GeecsDatabase.password)

    api_error.clear()
    device_ip, device_port = GeecsDatabase.find_device('U_ESP_JetXYZ')
    print(api_error)

    if device_ip:
        print('Device:\n\t' + device_ip + f', {device_port}')
    else:
        print('Device not found')
