import numpy as np
import astropy
import copy as copy
import os

from . import GenerateAverageModel as GAM
from . import SaveAverageOutput as SAO
from . import AverageModelPlot as AMP

"""
    This module averages the Barolo and FAT fits together to generate an average model.  It contains the routines:
    MakeAverageModelFolder --> This function makes a folder to store the average model files.
    AverageFits --> This function averages together the different fits.
    MakeAllMCGCubelets --> This function makes a set of MCG cubelets from the average model.
"""

def MakeAverageModelFolder(FittingOptions):
    """
         This function makes a folder to store the average model files.
    """
    os.makedirs(FittingOptions['AverageModelFolder'], exist_ok=True)

def AverageFits(ObjDict,FolderDict,AnalysisFncs,FittingOptions,FitFlags):
    """
        This function averages together the different fits.
    """
    #   First load in and analyze the sofia cubelet for comparisons later on (see CubeLoad/CubeAnalysis.py)
    CubeInfo=AnalysisFncs['CubeFnc']['SoFiA_CubeAnalysis'](ObjDict)
    #   Add in the object name from the cube file to the object dictionary
    #       Try to get the name of the galaxy from the cube header (but replace spaces with underscores)
    try:
        ObjDict['ObjName_From_Cube']=CubeInfo['CubeHeader']['OBJECT'].replace(' ','_')
    except:
        ObjDict['ObjName_From_Cube']=ObjDict['ObjName_Underscored']
    
    #   Get all the fits that will go into the averaging (see FitComparisons/CompareFits.py)
    FitParams=AnalysisFncs['FitCompFncs']['GetAllFits'](ObjDict,FolderDict,AnalysisFncs,FittingOptions,CubeInfo)
    #   Average the fits together into an average model
    AvgModel=GAM.MakeAvgModel(FitParams,FittingOptions,CubeInfo,AnalysisFncs,ObjDict)
    #   Once the average model is made, get the RA and DEC values and errors corresponding to the center point
    #       This needs to use the original, unconverted SoFiA cube header due to the fish-eye effect
    #           So load in the original cube
    OriCubeInfo=AnalysisFncs['CubeFnc']['BasicCubeAnalysis'](ObjDict['CubeFileName'])
    #   Now get the average model sky coordinates and errors
    AvgModel=GAM.GetAvgModelSkyCoordsAndErrs(AvgModel,OriCubeInfo,AnalysisFncs['AstroFncs'])
    #       Also, adjust the position angle as the RA-DEC grid is not strictly aligned with the pixel grid.
    AvgModel=GAM.GetGlobalPositionAngle(AvgModel,OriCubeInfo,AnalysisFncs['AstroFncs'])
    
    #   Correctly round the average model values
    AvgModel=GAM.RoundDict(AvgModel)
    #       Add the fit flags to the average model
    AvgModel['FitFlags']=FitFlags
    #   Save the Average Model
    AvgOutputDict=SAO.SaveAvgModel(AvgModel,ObjDict,FolderDict,FittingOptions,AnalysisFncs)
    #   Make and save all the MCG realizations that are needed -- and save the cubelet name to the ObjDict
    ObjDict=MakeAllMCGCubelets(ObjDict,AvgModel,AvgOutputDict,FittingOptions,AnalysisFncs)
    #   Compare the smoothed model cube to the smoothed data
    GAM.CompareAvgModelToFitData(AvgModel,ObjDict,AvgOutputDict,FittingOptions,AnalysisFncs['CubeFnc'])
    #   Copy the various fits
    SAO.CopyAllFits(ObjDict,FolderDict,FittingOptions,AvgOutputDict,AnalysisFncs)
    #   Copy the various data cubelets
    SAO.CopyAllData(ObjDict,AvgOutputDict,FittingOptions)
    #   Make a diagnostic plot
    AMP.MakeAvgModelPlot(ObjDict,FolderDict,CubeInfo,AvgModel,FittingOptions,AvgOutputDict,AnalysisFncs)
    #   Assign the object name, ID, and destination folder to the average model dictionary
    AvgModel['ID']=ObjDict['CatEntry']['id']
    AvgModel['NAME']=ObjDict['ObjName']
    AvgModel['FOLDER']=ObjDict['ObjName_Underscored']
    #   Return the average model for final catalogue construction
    return AvgModel
    

def MakeAllMCGCubelets(ObjDict,AvgModel,AvgOutputDict,FittingOptions,AnalysisFncs):
    """
        This function makes all the model cubelets necessary for the average model
    """
    #   When making the MCG cubelets, we need to convert back to Jy/beam km/s units
    AvgModel['SURFDENS_FACEON']*=AvgModel['SDCONV']
    #   Make an MCG cubelet at the full cube resolution
    SmoothSwitch=0
    MCGDict=GAM.MakeModelCube(AvgModel,AnalysisFncs['MCGFncs'],AnalysisFncs['CubeFnc'],ObjDict,SmoothSwitch)
    #   Save the cubelet and clean up the MCG directory
    ObjDict['FullResAvgModelCubeName']=SAO.MoveMCGModelRealization(ObjDict,MCGDict,AvgOutputDict,FittingOptions,SmoothSwitch)
    #   Make the MCG cubelet at the smoothed resolution
    SmoothSwitch=1
    MCGDict=GAM.MakeModelCube(AvgModel,AnalysisFncs['MCGFncs'],AnalysisFncs['CubeFnc'],ObjDict,SmoothSwitch)
    #   Save the cubelet and clean up the MCG directory
    ObjDict['SmoothedResAvgModelCubeName']=SAO.MoveMCGModelRealization(ObjDict,MCGDict,AvgOutputDict,FittingOptions,SmoothSwitch)
    #   Finally convert back to Msol/pc^2 units
    AvgModel['SURFDENS_FACEON']/=AvgModel['SDCONV']
    #   When converting back, we need to redo the rounding
    for i in range(len(AvgModel['SURFDENS_FACEON'])):
        AvgModel['SURFDENS_FACEON'][i]=round(AvgModel['SURFDENS_FACEON'][i],3)
    return ObjDict
