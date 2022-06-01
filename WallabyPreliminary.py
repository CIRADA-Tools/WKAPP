import numpy as np
import pandas as pd

import ReleaseConfigurationOptions as RCO
import SimpleIO as SIO
import PreprocessingRoutines as PPR
import BasicStatistics as BS

"""
    This script does a fair amount of the preprocessing necessary for the kinematic proto-pipeline.  It goes through every cube and converts them from frequency space to velocity space.  It also recenters the spatial coordinates of the headers in order to work with FAT and 3DBarolo.  However, this does introduce an error in the RA and Dec center (but not the pixel center) that must be accounted for afterwards.  Finally, each a 12 km/s velocity cube is made using 3DBarolo's smoothing function.
"""

def Main():
    #   Get the various file and folder names
    FolderDict=RCO.FF.BasicFolderAndFileLoc()
    #   Get the details about the excel file
    ExcelDetailsDict=RCO.FF.ExcelCatDetails()
    #   Load in the general WALLABY catalogue
    WallabyCat=SIO.SCIO.LoadXLS(FolderDict['CatalogueFile'],ExcelDetailsDict)
    #   Get the number of galaxies
    nGalaxies=np.shape(WallabyCat['id'])[0]
    #   Convert the ID column into integers
    WallabyCat['id']=WallabyCat['id'].astype('int')
    
    #   Set up the velocity folders
    PPR.VC.SetUpVelFolders(FolderDict)

    #for i in range(1,4):
    for i in range(1,nGalaxies+1):
        #   Get the basic object information
        ObjDict=RCO.FF.InitializeObjectDictionary(WallabyCat,i,FolderDict)
        #   Convert the cube header (both spatial and velocity)
        PPR.VC.CubeHeaderConvert(ObjDict,FolderDict,BS.AF.FuncDict())
    #   Smooth the cubelet
        PPR.BSS.BaroloSmoothing(ObjDict,FolderDict)

Main()
