import numpy as np
import pandas as pd

import ReleaseConfigurationOptions as RCO
import SimpleIO as SIO
import BasicStatistics as BS
import FATAnalysis as FP

"""
    This program is designed to set up the files needed to run FAT on all galaxies in a particular field with size > 2 beams or log(S/N)>1.25.  Unlike the Wallaby_BaroloAnalysis.py script, this one does not actually run FAT.
"""

def Main():

    #   Get the various file and folder names
    FolderDict=RCO.FF.BasicFolderAndFileLoc()
    #   Set the size limit
    SizeLimit=RCO.AO.SetSizeLimit()
    #   Set the FAT analysis options
    FATOptions=RCO.AO.FATAnalysisOptions(FolderDict)
    #   Get the details about the excel file
    ExcelDetailsDict=RCO.FF.ExcelCatDetails()
    #   Load in the general WALLABY catalogue
    WallabyCat=SIO.SCIO.LoadXLS(FolderDict['CatalogueFile'],ExcelDetailsDict)
    #   Get the number of galaxies
    nGalaxies=np.shape(WallabyCat['id'])[0]
    #   Convert the ID column into integers
    WallabyCat['id']=WallabyCat['id'].astype('int')
    
    #   Load in the Size and other measurements catalogue
    WallabyMeasurementCat=SIO.SCIO.LoadCSV(FolderDict['MeasurementCatalogue'],',')
    #   Get the FAT catalogue header
    FATCatStr=FP.FPS.FATHeaderStr()
    
    #   Make the FAT analysis folder
    FP.FPS.FATFolderCheck(FATOptions)
    #   Loop through all galaxies
    #       Keep a count of all galaxies in the FAT analysis catalogue
    Count=0
    #for i in range(1,5):
    for i in range(1,nGalaxies+1):
        #   Set up the specific galaxy dictionary object by the step
        ObjDict=RCO.FF.InitializeObjectDictionary(WallabyCat,i,FolderDict)
        #   Get the galaxy measurements from the measurement catalogue made by Wallaby_BasicAnalysis.py.  This is done by matching the ID.
        MeasureIndx=WallabyMeasurementCat.index[WallabyMeasurementCat['ID']==WallabyCat['id'][ObjDict['RowIndx']]][0]
        #   If the galaxy is large enough or has a high enough S/N, add it to the FAT analysis catalogue
        if WallabyMeasurementCat['SIZE'][MeasureIndx] >= SizeLimit or np.log10(WallabyMeasurementCat['SN_OBS'][MeasureIndx]) >= 1.25:
            #   The SetupFATObj function sets up a folder with a velocity cube for each galaxy in the FAT analysis folder and writes out a catalogue line needed for each galaxy in the actual FAT run.
            FATCatStr+= FP.FPS.SetupFATObj(ObjDict,FolderDict,FATOptions,Count)
        #   Increase the count for the FAT analysis
            Count=Count+1
            
    #   Write the FAT catalogue to a file
    FATCat=open(FATOptions['FATAnalysisCat'],'w')
    FATCat.write(FATCatStr)
    FATCat.close()
    
Main()
