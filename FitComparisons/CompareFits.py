import numpy as np
import astropy
import copy as copy
import os

from . import ComparisonPlot as CP
"""
    This module contains routines needed to compare a number of Barolo and FAT fits to a given galaxy with each other.  It contains the routines:
    FitComparisonFunctions --> This function makes a dictionary of other functions that may be useful in other modules.
    MakeFitComparisonsFolder --> This function makes a folder to hold the comparison plots.
    CompareFits --> This function does the actual fit comparisons.
    GetFits --> This function loads in all the different fits into an array of tilted ring dictionaries.
    GetBaroloFit --> This loads in a Barolo fit as a tilted ring dictionary
    GetFATFits --> This loads in a FAT fit as a pair of tilted ring dictionaries.
"""

def FitComparisonFunctions():
    """
        This function makes a dictionary of other functions for ease in passing to other modules.
    """
    CompareFncs={'GetAllFits':GetFits,'GetBaroloFit':GetBaroloFit,'GetFATFit':GetFATFits}
    return CompareFncs
    
def MakeFitComparisonsFolder(FolderDict):
    """
        This module makes a folder to holder the comparison plots
    """
    os.makedirs(FolderDict['FitsComparisionFolder'], exist_ok=True)

def CompareFits(ObjDict,FolderDict,AnalysisFncs,FittingOptions):
    """
        This function does the comparison between multiple Barolo and FAT fits and makes a diagnostic plot
    """
    #   First load in and analyze the sofia cubelet for comparisons later on.  The load-in routine is found in CubeLoad/CubeAnalysis.py
    CubeInfo=AnalysisFncs['CubeFnc']['SoFiA_CubeAnalysis'](ObjDict)
    #   Add in the object name from the cube file to the object dictionary
    #       Try to get the name of the galaxy from the cube header (but replace spaces with underscores)
    try:
        ObjDict['ObjName_From_Cube']=CubeInfo['CubeHeader']['OBJECT'].replace(' ','_')
        #   If this doesn't work, use the underscored name from the catalogue.
    except:
        ObjDict['ObjName_From_Cube']=ObjDict['ObjName_Underscored']
    #   Load all the fits into an array of Tilted-ring dictionaries
    FitParams=GetFits(ObjDict,FolderDict,AnalysisFncs,FittingOptions,CubeInfo)
    #   Make a comparison plot
    CP.MakeComparisonPlot(ObjDict,FolderDict,CubeInfo,FitParams,FittingOptions,AnalysisFncs)


def GetFits(ObjDict,FolderDict,AnalysisFncs,FittingOptions,CubeInfo):
    """
    This function loads in all the different fits into an array of tilted ring dictionaries.
    """
    #   Set the number of fit dictionaries that will be required
    FitParams=[None]*FittingOptions['nTotFits']
    #   Initialize the fit counter
    FitNum=0
    #   Loop through all Barolo Fits
    for i in range(FittingOptions['nBaroloFits']):
        #   Load in the Barolo fit
        FitParams[FitNum]=GetBaroloFit(ObjDict,FolderDict,AnalysisFncs,FittingOptions,CubeInfo,i)
        #   Increase the number of dictionaries by 1
        FitNum+=1
        
    #   Loop through all FAT fits
    for i in range(FittingOptions['nFATFits']):
        #   Load in the FAT fit
        FitParams[FitNum],FitParams[FitNum+1]=GetFATFits(ObjDict,FolderDict,AnalysisFncs,FittingOptions,CubeInfo,i)
        #   Increase the counter by 2 as each FAT fit has 2 dictionaries.
        FitNum+=2
    #   Return the array of dictionaries
    return FitParams


def GetBaroloFit(ObjDict,FolderDict,AnalysisFncs,FittingOptions,CubeInfo,FitStep):
    """
        This dictionary loads in a Barolo fit as a tilted ring dictionary
    """
    #   Load in the Barolo fit for the galaxy indicated by the ObjDict, using the particular fit that was run with the options set in the particular field (see ReleaseConfigurationOptions/).  The load-in functions are found in BaroloAnalyis/BaroloModelAnalysis.py
    BaroloFit=AnalysisFncs['BaroloFnc']['LoadBaroloModel'](ObjDict,FolderDict,FittingOptions['BaroloAnalysisFolders'][FitStep],FittingOptions['BaroloLabels'][FitStep])
        #   It is possible that the Barolo fit failed, and the load-in function returned the bad-fit dictionary.  If it didn't fail (indicated by the FitAchieved flag) there are a few more things to do
    if BaroloFit['FITAchieved']:
        #   Barolo gives the center in pixels, and some plots and calculations need this to be in RA and DEC, so get those values here.
        BaroloFit['RA'],BaroloFit['DEC']=AnalysisFncs['AstroFncs']['CalcRA_Dec_FromCube'](BaroloFit['XCENTER'],BaroloFit['YCENTER'],CubeInfo)
            #   If the fit was successful, compare an MCG realization of the parameters to the data cube
        BaroloFit=AnalysisFncs['ModelCompFncs']['CompareGeneralTiltedRingModel'](BaroloFit,CubeInfo,AnalysisFncs['MCGFncs'])
    #   Return the BaroloFit tilted-ring dictionary.
    return BaroloFit
    
def GetFATFits(ObjDict,FolderDict,AnalysisFncs,FittingOptions,CubeInfo,FitStep):
    """
        This dictionary loads in a FAT fit as a pair of tilted ring dictionaries
    """

    #   Start by making a 2 element empty list to store the pair of dictionaries
    FATFit=[None]*2
    #   Load in the FAT fit for the galaxy indicated by the ObjDict, using the particular fit that was run with the options set in the particular field (see ReleaseConfigurationOptions/).  The load-in functions are found in FATAnalyis/FATModelAnalysis.py
    FATFit[0],FATFit[1]=AnalysisFncs['FATFnc']['LoadFATModel'](ObjDict,FolderDict,FittingOptions['FATAnalysisFolders'][FitStep],FittingOptions['FATLabels'][FitStep])
    #   It is possible that the FAT fit failed, and the load-in function returned the bad-fit dictionary.  If it didn't fail (indicated by the FitAchieved flag) there are a few more things to do
    #       Loop through both dictionaries
    for j in range(2):
        if FATFit[j]['FITAchieved']:
            #   FAT returns the centroid in RA and DEC, which need to be converted to pixel X and Y values for some comparisons.
            FATFit[j]['XCENTER'],FATFit[j]['YCENTER']=AnalysisFncs['AstroFncs']['CalcCenter_FromCube'](FATFit[j]['RA'],FATFit[j]['DEC'],CubeInfo)
                #   If the fit is successful, compare an MCG realization of the parameters to the data cube
            FATFit[j]=AnalysisFncs['ModelCompFncs']['CompareGeneralTiltedRingModel'](FATFit[j],CubeInfo,AnalysisFncs['MCGFncs'])
    #   Return the pair of FAT Tilted Ring dictionaries.
    return FATFit
