import numpy as np
import pandas as pd


import ReleaseConfigurationOptions as RCO
import SimpleIO as SIO
import BaroloAnalysis as BA

"""
    This program is designed to run 3D Barolo on all galaxies in a particular field with size > 2 beams or log(S/N)>1.25
"""

def Main():

    #   Get the various file and folder names
    FolderDict=RCO.FF.BasicFolderAndFileLoc()
    #   Set the size limit for the analysis
    SizeLimit=RCO.AO.SetSizeLimit()
    #   Set the Barolo analysis options
    BaroloOptions=RCO.AO.BaroloAnalysisOptions(FolderDict)
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
    #   Set up a folder to contain the Barolo Analysis
    BA.BAS.BaroloAnalysisFolderPrep(BaroloOptions)

    #   Loop through all galaxies
    #for i in range(1,4):
    for i in range(1,nGalaxies+1):
        #   Set up the specific galaxy dictionary object by the step
        ObjDict=RCO.FF.InitializeObjectDictionary(WallabyCat,i,FolderDict)
        #   Get the galaxy measurements from the measurement catalogue made by Wallaby_BasicAnalysis.py.  This is done by matching the ID.
        MeasureIndx=WallabyMeasurementCat.index[WallabyMeasurementCat['ID']==WallabyCat['id'][ObjDict['RowIndx']]][0]
        #   If the galaxy is large enough or has a high enough S/N, run Barolo
        if WallabyMeasurementCat['SIZE'][MeasureIndx] >= SizeLimit or np.log10(WallabyMeasurementCat['SN_OBS'][MeasureIndx]) >= 1.25:
            BA.BAS.RunBaroloScript(ObjDict,FolderDict,BaroloOptions)
            

Main()
