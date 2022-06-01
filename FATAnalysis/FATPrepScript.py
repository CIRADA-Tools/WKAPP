import numpy as np
import astropy
import copy as copy
import os

"""
    This module contains routines for running setting up (but not running) the FAT analysis on the WALLABY  galaxies.  It contains the routines:
    FATHeaderStr --> This routine sets up the header for a FAT catalogue file
    FATFolderCheck --> This file makes a folder to organize the files needed for the FAT analysis
    SetupFATObj --> This routine sets up a folder for a specific galaxy for the FAT analysis as well as writing out the string needed for the FAT catalogue file.
"""

def FATHeaderStr():
    """
        This routine writes the header for a FAT catalogue.
    """
    HeaderStr="number|Distance|Directoryname|Cubename\n"
    return HeaderStr
    
def FATFolderCheck(FATOptions):
    """
        This routine makes a folder for the full FAT analysis
    """
    os.makedirs(FATOptions['FATAnalysisFolder'], exist_ok=True)

def SetupFATObj(ObjDict,FolderDict,FATOptions,Count):
    """
        This routine sets up the folder structure for a specific galaxy in a FAT run and writes out the string required in the FAT catalogue file.
    """

    #   Get the velocity cube name.  Use the smoothing switch set in the field configuration routines (see /ReleaseConfigurationOptions/ for the routines where this is set) to decide whether to use the 4 km/s or 12 km/s cubes.
    if FATOptions['UseSmoothedCubes']:
        TargCube=ObjDict['SmoothedVelocityCubeFileName']
    else:
        TargCube=ObjDict['VelocityCubeFileName']
    #   It's necessary to separate the specific cube name from the larger folder structure.
    CubeName=TargCube.split("/")[-1]
    #   Name the analysis folder that will live inside the larger FAT directory
    AnalysisFolder=FolderDict['BaseName']+"_"+ObjDict['NumStr']
    #   Get the full path for the galaxy analysis
    ObjFolder=FATOptions['FATAnalysisFolder']+AnalysisFolder+"/"
    #   Make the Analysis folder
    os.makedirs(ObjFolder,exist_ok=True)
    #   Copy the velocity cube to the folder
    os.system("cp "+ TargCube+ " " + ObjFolder+".")
    #   Make out the catalogue string needed for running FAT
    #       The catalogue needs a name without a .fits suffix.
    pyFatCubeName=CubeName.split(".")[0]
    #   The catalogue file follows the format Count|Distance|Folder|CubeName
    #       Using -1 as the distance tells FAT to calculate it internally.
    CatStr=str(Count)+"|-1|"+AnalysisFolder+"|"+pyFatCubeName +"\n"
    return CatStr

