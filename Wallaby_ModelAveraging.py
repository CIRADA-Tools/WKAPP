import numpy as np
import pandas as pd

from collections import Counter

import ReleaseConfigurationOptions as RCO
import SimpleIO as SIO
import BaroloAnalysis as BA
import FATAnalysis as FA
import BasicStatistics as BS
import DiagnosticPlots as DP
import MCGModel as MM
import CompareModel as CM
import CubeLoad as CL

import FitComparisons as FC
import ModelAveraging as MA
"""
    This module generates a pilot phase kinematic model by averaging together a set of FAT and 3DBarolo fits.
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
    ComparisonFncs=FC.CF.FitComparisonFunctions()
    AnalysisFncs={'CubeFnc':CubeFnc,'BaroloFnc':BaroloFnc,'FATFnc':FATFnc,'AstroFncs':AstroFncs,'PlotFncs':PlotFncs,'MCGFncs':MCGFncs,'ModelCompFncs':ModelCompFncs,'FitCompFncs':ComparisonFncs}
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
    #   Load in the Fit success catalogue
    FitSuccessCat=SIO.SCIO.LoadCSV(FolderDict['FitSuccessCatalogue'],',')
    #   Set up a folder to contain the fit comparison plots
    MA.AF.MakeAverageModelFolder(FittingOptions)
    #   Figure out how many successful models there are from the fit success catalogue (made by hand after running Wallaby_BasicFitComparison)
    nUsuable=Counter(FitSuccessCat['USABLE_MODEL'])[1]
    #   Set up an array for the average models
    AvgModels=[None]*nUsuable
    #   Start a counter for the number of kinematic models generated
    Tot=0
    #   Loop through all galaxies
    #for i in range(1,4):
    for i in range(1,nGalaxies+1):
        #   Set up the specific galaxy dictionary object by the step
        ObjDict=RCO.FF.InitializeObjectDictionary(WallabyCat,i,FolderDict)
        #   Print out the current status
        print("Getting Fits for step ", i, " which is ", ObjDict['ObjName'])
        #   Get the galaxy measurements from the measurement catalogue made by Wallaby_BasicAnalysis.py.  This is done by matching the ID.
        MeasureIndx=WallabyMeasurementCat.index[WallabyMeasurementCat['ID']==WallabyCat['id'][ObjDict['RowIndx']]][0]
        #   Get the Barolo and FAT fit success tags from the fit success catalogue by matching the ID
        FitIndx=FitSuccessCat.index[FitSuccessCat['ID']==WallabyCat['id'][ObjDict['RowIndx']]][0]
        #   Copy the specifc row of the measurements into the galaxy dictionary
        ObjDict['MeasureCat']=WallabyMeasurementCat.loc[WallabyMeasurementCat['ID']==WallabyCat['id'][ObjDict['RowIndx']]]
        #   If the fits were judged by eye to be a success, generate an average model
        if FitSuccessCat['USABLE_MODEL'][FitIndx] ==1 :
            #   Print a note to average a model
            print("Average this model", ObjDict['ObjName'],FitSuccessCat['USABLE_MODEL'][FitIndx])
            #   Make the average model and store the parameters into the array of successful tilted ring models.
            AvgModels[Tot]=MA.AF.AverageFits(ObjDict,FolderDict,AnalysisFncs,FittingOptions,FitSuccessCat['FLAGS'][FitIndx])
            #   Keep track of the total number of successful models
            Tot+=1
    print("Total successful models", Tot)
    #   Save the average models to a csv file
    MA.SAO.WriteOutputTable(AvgModels,FittingOptions)

Main()
