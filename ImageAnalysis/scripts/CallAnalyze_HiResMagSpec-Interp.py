# -*- coding: utf-8 -*-

"""
Very simple front end that showcases how the HiResMagSpec image analyzer is called.

@auther Chris
"""

# Imports
import sys
sys.path.insert(0, "../")
import modules.MagSpecAnalysis as MagSpecAnalysis
import modules.DirectoryModules as DirectoryFunc

"""
Generate the filepath to the saved data here.  Including the scan and shot number.
In the future, this step will be embedded in some other process.
"""

sampleCase = 2
if sampleCase == 1:
    data_day = 29
    data_month = 6
    data_year = 2023

    scan_number = 23
    shot_number = 1

elif sampleCase == 2:
    data_day = 13
    data_month = 7
    data_year = 2023

    scan_number = 20
    shot_number = 5

else:
    print("Pick a valid sample case!")
    data_day = 0; data_month = 0; data_year = 0

superpath = DirectoryFunc.CompileDailyPath(data_day, data_month, data_year)
image_name = "U_HiResMagCam-interp"

"""
From the path, load the image, tdms, and interpSpec files
LoadImage by default includes background reduction and charge normalization.  I turn it off here.
"""

raw_image = MagSpecAnalysis.LoadImage(superpath, scan_number, shot_number, image_name, doThreshold = False, doNormalize = False)
tdms_filepath = DirectoryFunc.CompileTDMSFilepath(superpath, scan_number)
interpSpec_filepath = DirectoryFunc.CompileFileLocation(superpath, scan_number, shot_number,
                                                        imagename = 'U_HiResMagCam-interpSpec', suffix=".txt")

"""
Below is some optional code to test that I can skip right to the plots and not worry about the dictionary.
"""

#MagSpecAnalysis.PlotBeamDistribution(raw_image, tdms_filepath, interpSpec_filepath, shot_number, doThreshold = True, doNormalize = True)
#sys.exit()

"""
Call the overhead "analyze" function with all of this, and return a dictionary of parameters

Note:  this function first performs background reduction and then charge normalization.
"""

MagSpecDict = MagSpecAnalysis.AnalyzeImage_Interp(raw_image, tdms_filepath, interpSpec_filepath, shot_number)
print(MagSpecDict)

"""
For my own benefit, call a separate function to make the nice plots
"""

plotInfo = DirectoryFunc.CompilePlotInfo(data_day, data_month, data_year, scan_number, shot_number, image_name)
MagSpecAnalysis.PlotEnergyProjection_Interp(raw_image, tdms_filepath, interpSpec_filepath, shot_number, plotInfo = plotInfo, analyzeDict = MagSpecDict, doThreshold = True, doNormalize = True)
MagSpecAnalysis.PlotSliceStatistics_Interp(raw_image, tdms_filepath, interpSpec_filepath, shot_number, plotInfo = plotInfo, analyzeDict = MagSpecDict, doThreshold = True, doNormalize = True)
MagSpecAnalysis.PlotBeamDistribution_Interp(raw_image, tdms_filepath, interpSpec_filepath, shot_number, plotInfo = plotInfo, analyzeDict = MagSpecDict, doThreshold = True, doNormalize = True)
