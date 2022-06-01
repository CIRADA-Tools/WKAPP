"""
    This module contains the definitions needed to analyze the NGC 4636 field sources.  It contains the routines:
    FolderAndFileLoc --> This sets the locations of folders and files for the Hydra field
    ExcelCatDetails --> This sets the definitions of the SoFiA excel catalogue file
    FATAnalysisOptions --> This sets the definitions for the FAT analysis
    BaroloAnalysisOptions --> This sets the definitions for the Barolo analysis
    FitComparisonOptions --> This sets the definitions for comparing and combining the FAT and Barolo analyses
"""
def FolderAndFileLoc(ParentFolder):
    """
        This function gives the names of files and folders used in the various analysis scripts
    """
    #   Set the name of the source
    SourceName="WALLABY NGC4636 DR1"
    #   Set the name of the folder where the analysis will be done
    DataReleaseFolder=ParentFolder+"NGC4636DR1Analysis/"
    #   Set the name of the folder containing the sources
    SourceFolder=DataReleaseFolder+SourceName+"/"
    #   Set the base name of files in the release
    BaseName="WALLABY_PS_NGC4636_DR1"
    #   Set the name of the SoFiA catalogue file
    CatalogueFile=SourceFolder+BaseName+"_catalog.xls"
    #   Set the name of the folder containing the SoFiA outputs
    SofiADataFolder=SourceFolder+BaseName+"_source_products/"
    #   Set the name of the folder that will contain the velocity converted cubes
    VelFolder=DataReleaseFolder+"NGC4636_Vel_Cubes_Jan2022/"
    #   Set the name of the folder that will contain the smoothed velocity cubelets
    SmoothedVelFolder=VelFolder+"chan_12kms/"
    #   Set the location of the Barolo executable
    BaroloPath="/Users/nate/Dropbox/Bbarolo-1.6/BBarolo_1.6_Precompiled"
    #   Set the name of the folder that will contain all the Barolo analyses -- this helps keep things organized if Barolo is run more than once
    MainBaroloFolder=DataReleaseFolder+"Barolo_Analysis/"
    #   Set the name of the Folder that will contain all FAT analyses -- this helps keep things organized if Barolo is run more than once
    MainFATFolder=DataReleaseFolder+"FAT_Analysis/"
    #   Name the csv file that will contain catalogue measurements
    MeasurementCatalogue=DataReleaseFolder+BaseName+"_CatalogueMeasurements.csv"
    #   Name the location for the fit comparison plots
    FitsComparisionFolder=DataReleaseFolder+"FitComparions_Mar2022/"
    #   Set the name of the file that will contain the fit success flags.  This file, unlike most of the others, must be made by hand
    FitSuccessCatalogue=DataReleaseFolder+BaseName+"_FitSuccess_Flags.csv"
    #   Place all the definitions into a dictionary
    FolderDict=locals()
    return FolderDict

def ExcelCatDetails():
    """
        This function defines the details of the excel file containing the SoFiA catalogue.
    """
    #   Name the sheet
    SheetName="NGC4636"
    #   Give the length of the header
    HeaderLen=7
    #   Give the column that will be used for defining the dataframe
    StripCol=[1,'id']
    #   Place everything into a dictionary
    ExcelDetailsDict=locals()
    return ExcelDetailsDict
    
def FATAnalysisOptions(FolderDict):
    """
        This function contains definitions for the FAT analysis
    """
    #   Name the folder to hold the set of FAT analysis results
    FATAnalysisFolder=FolderDict['MainFATFolder']+"Jan2022_LowSize/"
    #   Name the catalogue file that FAT will need
    FATAnalysisCat=FATAnalysisFolder+FolderDict['BaseName']+"_FatCatalogue_LowSize.txt"
    #   Select whether to use smoothed or unsmoothed cubelets.
    UseSmoothedCubes=True
    #   Place the definitions into a dictionary
    FATOptions={'FATAnalysisFolder':FATAnalysisFolder,'FATAnalysisCat':FATAnalysisCat,'UseSmoothedCubes':UseSmoothedCubes}
    return FATOptions

def BaroloAnalysisOptions(FolderDict):
    """
        This function contains definitions for the Barolo analysis
    """
    #   Name the specific folder to hold the results of the Barolo run
    SpecificAnalysisFolder="Mar_2022_Analysis/"
    #   Select whether to use smoothed or unsmoothed cubelets.
    UseSmoothedCubes=True
    #   Select whether to use SoFiA masks or the inbuilt masking software
    UseSoFiAMasks=False
    #   Select whether to use SoFiA masks or the inbuilt masking software
    nFitsPerGalaxy=1
    #   Make an empty list to for the general names of the Barolo runs and the specific folders for each run
    FitNames=[None]*nFitsPerGalaxy
    BaroloAnalysisFolders=[None]*nFitsPerGalaxy
    #   Loop through all Barolo fits
    for i in range(nFitsPerGalaxy):
        if i == 0:
            #   Set the specific name of the Barolo analysis
            FitNames[i]="NoIni"
            #   Set the folder name for the Barolo run
            BaroloAnalysisFolders[i]=FolderDict['MainBaroloFolder']+SpecificAnalysisFolder+FitNames[i]+"/"
    #   Store all definitions into a dictionary
    BaroloOptions={'BaroloAnalysisFolders':BaroloAnalysisFolders,'UseSmoothedCubes':UseSmoothedCubes,'FitNames':FitNames,'nFits':nFitsPerGalaxy,'SpecificAnalysisFolder':SpecificAnalysisFolder,'UseSoFiAMasks':UseSoFiAMasks}
    return BaroloOptions

def FitComparisonOptions(FolderDict):
    """
        This function sets the definitions for the Barolo-FAT comparisons and for generating an average model from different fits
    """
    #   First set the version of the models
    ModelVersion=1
    #   Name the average model release
    AvgModelReleaseBaseName="Wallaby_NGC4636_DR1_KinematicModels_v"+str(ModelVersion)
    #   Set the average model release folder
    AverageModelFolder=FolderDict['DataReleaseFolder']+AvgModelReleaseBaseName+"/"
    #   Set the number of individual FAT and Barolo fits to be used in the averaging
    nBaroloFits=1
    nFATFits=1
    #   Set the total number of fits/TR dictionaries involved
    nTotFits=nBaroloFits+2*nFATFits
    #   Set up the Barolo fitting options
    #       First make empty lists for the folder names, label, fit names, and output names
    BaroloAnalysisFolders=[None]*nBaroloFits
    BaroloLabels=[None]*nBaroloFits
    BaroloOutFitName=[None]*nBaroloFits      #This is for naming the fit folders in the Avg model
    BaroloFitNames=[None]*nBaroloFits
        #   Loop through all the Barolo fits
    for i in range(nBaroloFits):
        if i ==0:
            #   Set the specific analysis folder containing the Barolo fit
            SpecificAnalysisFolder="Mar_2022_Analysis/"
            #   Set the name of the fit
            BaroloFitNames[i]="NoIni"
            #   Set the name of the folder containing the fits
            BaroloAnalysisFolders[i]=FolderDict['MainBaroloFolder']+SpecificAnalysisFolder+"/NoIni/"
            #   Set the label for the Barolo run in the comparison plots
            BaroloLabels[i]="Barolo - No Ini"
            #   Set the name of Barolo output names
            BaroloOutFitName[i]="Barolo-No_Ini"
    #   Do the same for FAT
    #       Again, make the empty lists for the folders, labels, and base configuration files
    FATAnalysisFolders=[None]*nFATFits
    FATLabels=[None]*nFATFits
    FATConfigFileName=[None]*nFATFits
    #   Loop through all the FAT fits
    for i in range(nFATFits):
        if i ==0:
            #   For each fit, give 2 labels for the approaching and receding dictionaries
            FATLabels[i]=["FAT D1","FAT"]
            #   Set the name of the FAT folder containing the results
            FATAnalysisFolders[i]=FolderDict['MainFATFolder']+"Jan2022_Sample/"
            #   Set the name of the base configuration file
            FATConfigFileName[i]="FAT_NGC4636_12kms.config"
    #   Add in the provenance key values
    SBID="10809,10812,10736"
    SRCVER="SoFiA 2.3.1"
    SRCTR="NGC4636 TR1"
    KINVER="WKAPP v1"
    KINTR="NGC4636 Kin TR1"
    ProvenanceKeyVals={'SBID':SBID,'SRCVER':SRCVER,'SRCTR':SRCTR,'KINVER':KINVER,'KINTR':KINTR}
    #   Store all the names and configuration switches into a dictionary
    FittingOptions={'nTotFits':nTotFits,'nBaroloFits':nBaroloFits,'nFATFits':nFATFits,'BaroloAnalysisFolders':BaroloAnalysisFolders,'FATAnalysisFolders':FATAnalysisFolders,'BaroloLabels':BaroloLabels,'BaroloOutFitName':BaroloOutFitName,'BaroloFitNames':BaroloFitNames,'FATLabels':FATLabels,'FATConfigFileName':FATConfigFileName,'ModelVersion':ModelVersion,'AverageModelFolder':AverageModelFolder,'AvgModelReleaseBaseName':AvgModelReleaseBaseName,'ProvenanceVals':ProvenanceKeyVals}
    return FittingOptions

    return FittingOptions
