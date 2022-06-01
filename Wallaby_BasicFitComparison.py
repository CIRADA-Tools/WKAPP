import numpy as np
import pandas as pd

import ReleaseConfigurationOptions as RCO
import SimpleIO as SIO
import BaroloAnalysis as BA
import FATAnalysis as FA

import FitComparisons as FC
import BasicStatistics as BS
import DiagnosticPlots as DP
import MCGModel as MM
import CompareModel as CM
import CubeLoad as CL

"""
    This module compares a set of Barolo and FAT fits to each other and makes a set of diagnostic plots that are then used to determine which fits are sufficiently successful to be averaged together.
"""


def Main():
    #   Load in the functions that will be needed
    CubeFnc=CL.CA.CubeAnalysisFuncDict()
    BaroloFnc=BA.BMA.BaroloModelFnc()
    FATFnc=FA.FMA.FATModelFnc()
    AstroFncs=BS.AF.FuncDict()
    PlotFncs=DP.DPF.PlotFnc()
    MCGFncs=MM.MMG.MCGModelFnc()
    ModelCompFncs=CM.CMTC.ComparisonFncs()
    AnalysisFncs={'CubeFnc':CubeFnc,'BaroloFnc':BaroloFnc,'FATFnc':FATFnc,'AstroFncs':AstroFncs,'PlotFncs':PlotFncs,'MCGFncs':MCGFncs,'ModelCompFncs':ModelCompFncs}
    #   Get the various file and folder names
    FolderDict=RCO.FF.BasicFolderAndFileLoc()
    #   Set the size limit for the analysis
    SizeLimit=RCO.AO.SetSizeLimit()
    #   Set the fit comparions/analysis options
    FittingOptions=RCO.AO.FitComparisonOptions(FolderDict)
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
    #   Set up a folder to contain the fit comparison plots
    FC.CF.MakeFitComparisonsFolder(FolderDict)
    #   Loop through all the galaxies
    #for i in range(1,4):
    for i in range(1,nGalaxies+1):
        #   Set up the specific galaxy dictionary object by the step
        ObjDict=RCO.FF.InitializeObjectDictionary(WallabyCat,i,FolderDict)
        print("Getting Fits for step ", i, " which is ", ObjDict['ObjName'])
        #   Get the galaxy measurements from the measurement catalogue made by Wallaby_BasicAnalysis.py.  This is done by matching the ID.
        MeasureIndx=WallabyMeasurementCat.index[WallabyMeasurementCat['ID']==WallabyCat['id'][ObjDict['RowIndx']]][0]
        #   Copy the specifc row of the measurements into the galaxy dictionary
        ObjDict['MeasureCat']=WallabyMeasurementCat.loc[WallabyMeasurementCat['ID']==WallabyCat['id'][ObjDict['RowIndx']]]
        #   If the galaxy is large enough or has a high enough S/N, make a comparison plot.
        if WallabyMeasurementCat['SIZE'][MeasureIndx] >= SizeLimit or np.log10(WallabyMeasurementCat['SN_OBS'][MeasureIndx]) >= 1.25:
            FC.CF.CompareFits(ObjDict,FolderDict,AnalysisFncs,FittingOptions)
            


Main()
