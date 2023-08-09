"""
Mon 8-7-2023

Module for containing functions to both calculate the energy axis calibration and retrieve the previously found fit
parameters.  Works for both the HiResMagSpec and the BCaveMagSpec

@Chris
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy import optimize
from scipy.optimize import curve_fit
#import CedossMathTools as MathTools

def return_default_energy_axis(pixel_axis):
    """
    This is calibrated using the 825 mT settings to get a range of 89.368 to 114.861 MeV across 1287 pixels
    """
    energy_axis = np.zeros(len(pixel_axis))
    magField = '825mT'
    if magField == '800mT':
        p0 = 8.66599527e+01
        p1 = 1.78007126e-02
        p2 = 1.10546749e-06
    if magField == '825mT':
        p0 = 8.93681013e+01
        p1 = 1.83568540e-02
        p2 = 1.14012869e-06
    energy_axis = p0 + p1*pixel_axis + p2*np.power(pixel_axis, 2)
    return energy_axis


def read_double_array(file_path):
    try:
        double_array = np.loadtxt(file_path, delimiter='\t')
        return double_array
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return np.array([])
    except Exception as e:
        print(f"An error occurred: {e}")
        return np.array([])

def polyfit_function(axis, fit):
    func = np.polyval(fit, axis)
    label = 'Order ' + str(len(fit)-1) + ': '
    for i in range(len(fit)):
        power = len(fit) - i - 1
        label += "{:.1e}".format(fit[i]) + '*p^' + str(power)
        if i != (len(fit) - 1):
            label += ' + '
    return func, label

def fitfunction(axis, data, function, guess=[0,0,0]):
    p0 = guess  # Adjust initial parameter guess based on your data
    popt, pcov = curve_fit(function, axis, data, p0=p0,
                           bounds=([-np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf]))
    return popt

def exponential(x, *p):
    a, b, c = p
    return np.exp(a * (1e-6*x) - b) + c


if __name__ == "__main__":
    mode = 'triple'
    if mode == 'triple':
        superpath = 'Z:/data/Undulator/Y2023/08-Aug/23_0802/auxiliary data/'
        magSetting = '825'
        cam1_path = superpath + 'magCam1/' + 'energy vs pixel MagCam1 ' + magSetting + ' mT.tsv'
        cam2_path = superpath + 'magCam2/' + 'energy vs pixel MagCam2 ' + magSetting + ' mT.tsv'
        cam3_path = superpath + 'magCam3/' + 'energy vs pixel MagCam3 ' + magSetting + ' mT.tsv'

        separation_pixels = 50

        cam1_array = read_double_array(cam1_path)
        cam1_axis = np.arange(0, len(cam1_array))
        cam2_array = read_double_array(cam2_path)
        cam2_axis = np.arange(0, len(cam2_array))+cam1_axis[-1]
        cam3_array = read_double_array(cam3_path)
        cam3_axis = np.arange(0, len(cam3_array))+cam2_axis[-1]+separation_pixels

        cam1_axis = cam1_axis
        cam2_axis = cam2_axis
        cam3_axis = cam3_axis

        section123_array = np.append(np.append(cam1_array, cam2_array), cam3_array)
        section123_axis = np.append(np.append(cam1_axis, cam2_axis), cam3_axis)
        section123_plotaxis = np.linspace(0,section123_axis[-1],100)
        section123_fit = np.polyfit(section123_axis, section123_array, 12)
        section123_func, section123_label = polyfit_function(section123_plotaxis, section123_fit)

        section12_array = np.append(cam1_array, cam2_array)
        section12_axis = np.append(cam1_axis, cam2_axis)
        section12_fit = np.polyfit(section12_axis, section12_array, 2)
        section12_func, section12_label = polyfit_function(section123_plotaxis, section12_fit)

        section3_fit = np.polyfit(cam3_axis, cam3_array, 7)
        section3_func, section3_label = polyfit_function(section123_plotaxis, section3_fit)

        #section3_efit = fitfunction(cam3_axis, cam3_array, exponential, guess=[1, 0.01, 0])
        #section3_efunc = exponential(section123_plotaxis, *section3_efit)

        #section3_ifit = MathTools.FitDataSomething(cam3_array, cam3_axis, MathTools.GaussianOffset, guess=[2000, 0.3, 1, 200])
        #section3_ifunc = MathTools.GaussianOffset(section3_ifit, section123_plotaxis)
        #print(section3_ifit)

        plt.plot(cam1_axis, cam1_array, label='cam1')
        plt.plot(cam2_axis, cam2_array, label='cam2')
        plt.plot(cam3_axis, cam3_array, label='cam3')
        plt.plot(section123_plotaxis, section12_func, ls = 'dotted', label=section12_label)
        #plt.plot(section123_plotaxis, section123_func, ls = 'dotted', label=section123_label)
        plt.plot(section123_plotaxis, section3_func, ls = 'dotted', label=section3_label)
        #plt.plot(section123_plotaxis, section3_efunc, ls ='dotted', label="exponential")
        #plt.plot(section123_plotaxis, section3_ifunc, ls ='dashed', label="gaussian")
        plt.ylim([min(cam1_array)*0.7,max(cam3_array)*1.1])
        plt.xlim([section123_plotaxis[0], section123_plotaxis[-1]])
        plt.legend()
        plt.show()


    if mode == 'single':
        superpath = 'Z:/data/Undulator/Y2023/08-Aug/23_0802/auxiliary data/HiResMagCam/'
        filename = 'energy vs pixel HiResMagCam 825 mT.tsv'
        #filename = 'energy vs pixel HiResMagCam 800 mT.tsv'
        file_path = superpath + filename
        double_array = read_double_array(file_path)
        if double_array.size > 0:
            print("Array of doubles:", double_array)
            axis = np.arange(0, len(double_array))
            fit = np.polyfit(axis, double_array, 2)
            print("Fit: (0th order last)")
            print(" ", fit)

            func, label = polyfit_function(axis, fit)

            plt.plot(axis, double_array, label = filename)
            plt.plot(axis, func, ls='--', label = label)
            plt.legend()
            plt.show()

