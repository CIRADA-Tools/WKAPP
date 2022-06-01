import numpy as np
import astropy
import os
import pandas as pd
from datetime import date

import astropy.units as u
from astropy.io import fits

"""
    This module contains routines for saving the outputs from the average model generation.  It contains the routines:
    AvgOutputFncDict --> This places useful routines into a dictionary for ease of access
    NameOutputs --> This function names the various output files/directories
    SaveAvgModel --> This function saves all data products for the average model.
    MakeOutputFolderStructure --> This function makes the output folder structure.
    WriteAvgModel --> This function writes the average model text file
    WriteRCFile --> This function writes the rotation curve to a fits file
    WriteSDFile --> This function writes the surface density to a fits file
    WriteGeoFile --> This function writes the model geometry to a fits file
    
    CopyAllFits --> This function copies the FAT and Barolo fit information into the results folder
    CopyFATFits --> This function copies in the FAT results
    CopyBaroloFits --> This function copies in the Barolo results
    CopyAllData --> This function copies in the data cubes that were used in the fitting process.
    MoveMCGModelRealization --> This copies the MCG realizations of the average model into the appropriate folders.
    
    AdjustHeaderPosBack --> This function is used to adjust the datacube header spatial reference point back to the SoFiA detection to avoid fisheye effects
    AdjustHeaderVelBack --> This function is used to adjust the datacube referece velocity to match the smoothing results
    
    WriteOutputTable --> This function writes all the models to a single csv file
    MakePANDASDataFrame --> This function turns the list of TR dictionaries into a PANDAS data frame
    
"""

def AvgOutputFncDict():
    """
       This function contains a dictionary of useful functions that can be utilized elsewhere
    """
    AvgOutFncDict={'NameOutputs':NameOutputs}
    return AvgOutFncDict
    
def NameOutputs(ObjDict,FittingOptions):
    """
        This function defines the names of various output files/folders
    """
    #   Name the main output folder
    OutputDirName=FittingOptions['AverageModelFolder']+ObjDict['ObjName_Underscored']+"/"
    #   Name the folder that will contain the individual fits
    FitsDir=OutputDirName+"FitsUsedForAveraging"
    #   Name the folder that will contain the full resolution cubes
    FullResDir=OutputDirName+"FullResolution"
    #   Name all the FAT Fit folders
    FATDir=[None]*FittingOptions['nFATFits']
    for i in range(FittingOptions['nFATFits']):
        FATDir[i]=FitsDir+"/"+FittingOptions['FATLabels'][i][1]
    #   Name all the Barolo Fit folders
    BaroloDir=[None]*FittingOptions['nBaroloFits']
    for i in range(FittingOptions['nBaroloFits']):
        BaroloDir[i]=FitsDir+"/"+FittingOptions['BaroloOutFitName'][i]
    #   Set the base name for output files
    KinTR_Underscore=FittingOptions['ProvenanceVals']['KINTR'].split()
    KinTR_Underscore=KinTR_Underscore[0]+"_"+KinTR_Underscore[1]+"_"+KinTR_Underscore[2]
    OutputBaseName=ObjDict['ObjName_Underscored']+"_"+KinTR_Underscore
    #   Store everything in a dictionary and return it
    OutputNameDict={'OutputDirName':OutputDirName,'FitsDir':FitsDir,'FullResDir':FullResDir,'FATDir':FATDir,'BaroloDir':BaroloDir,'OutputBaseName':OutputBaseName}
    return OutputNameDict

def SaveAvgModel(AvgModel,ObjDict,FolderDict,FittingOptions,AnalysisFncs):
    """
        This function saves all the data products that go into a WALLABY pilot phase I kinematic model
    """
    #   Name all the output folders
    AvgOutputDict=NameOutputs(ObjDict,FittingOptions)
    #   Make the main output folder
    MakeOutputFolderStructure(AvgOutputDict,FittingOptions)
    #   Write up the average model parameter file
    WriteAvgModel(AvgModel,ObjDict,AvgOutputDict,FittingOptions,FolderDict)
        #   Write a fits file for the rotation curve
    WriteRCFile(AvgModel,ObjDict,AvgOutputDict,FittingOptions,FolderDict)
        #   Write a fits file for the surface density
    WriteSDFile(AvgModel,ObjDict,AvgOutputDict,FittingOptions,FolderDict)
        #   Write a fits file for the geometric parameters
    WriteGeoFile(AvgModel,ObjDict,AvgOutputDict,FittingOptions,FolderDict)
    #   Return the output dictionary for use elsewhere
    return AvgOutputDict
    
def MakeOutputFolderStructure(OutputDict,FittingOptions):
    """
        This function makes the folders needed for a particular galaxy's kinematic model
    """
    #   Make the main output directory
    os.makedirs(OutputDict['OutputDirName'], exist_ok=True)
    #   Make the fits comparison folder
    os.makedirs(OutputDict['FitsDir'], exist_ok=True)
    #   Make the full resolution folder
    os.makedirs(OutputDict['FullResDir'], exist_ok=True)
    return
    
def CopyAllFits(ObjDict,FolderDict,FittingOptions,AvgOutputDict,AnalysisFncs):
    """
        This function copies in the both the FAT and Barolo results to the appropriate folders.
    """
      #   Save the FAT fits
    CopyFATFits(ObjDict,FolderDict,FittingOptions,AvgOutputDict,AnalysisFncs['FATFnc'])
    #   Save the Barolo fits
    CopyBaroloFits(ObjDict,FolderDict,FittingOptions,AvgOutputDict,AnalysisFncs['BaroloFnc'])
    
def CopyAllData(ObjDict,AvgOutputDict,FittingOptions):
    """
        This function copies the data used in the fitting process into the appropriate folders.
    """
    #   First set the target name for the smoothed cube
    SmoothedCubeDestName=AvgOutputDict['OutputDirName']+AvgOutputDict['OutputBaseName']+"_ProcessedData.fits"
    #   Next set the target name for the full resolution cubelet
    FullResCubeDestName=AvgOutputDict['FullResDir']+"/"+AvgOutputDict['OutputBaseName']+"_FullResProcessedData.fits"
    #   Copy smoothed cube
    os.system("cp "+ObjDict['SmoothedVelocityCubeFileName']+" "+ SmoothedCubeDestName)
    #   Copy full resolution cube
    os.system("cp "+ObjDict['VelocityCubeFileName']+" "+ FullResCubeDestName)
    #   In order to avoid fish-eye issues reset the spatial header to the reference points used in the original SoFiA cubelets
    AdjustHeaderPosBack(SmoothedCubeDestName,ObjDict['CubeFileName'])
    AdjustHeaderPosBack(FullResCubeDestName,ObjDict['CubeFileName'])
    #   Add the provenance key words to the files
    AddProvenanceKeyWords(SmoothedCubeDestName,FittingOptions['ProvenanceVals'],0,ObjDict)
    AddProvenanceKeyWords(FullResCubeDestName,FittingOptions['ProvenanceVals'],0,ObjDict)
def MoveMCGModelRealization(ObjDict,MCGDict,AvgOutputDict,FittingOptions,SmoothSwitch):
    """
        This function copies the MCG realizations of the kinematic TR models into the appropriate folders
    """
    #   Set the destination folder and cube name
    #       Smoothswitch=0 is full resolution
    if SmoothSwitch == 0:
        DestFolder=AvgOutputDict['FullResDir']
        DestCubeName=AvgOutputDict['OutputBaseName']+"_FullResModelCube.fits"
    #       Smoothswitch=1 is at the smoothed cubelet resolution
    if SmoothSwitch == 1:
        DestFolder=AvgOutputDict['OutputDirName']
        DestCubeName=AvgOutputDict['OutputBaseName']+"_ModelCube.fits"
    #   Set the name of the cube
    FullDestName=DestFolder+"/"+DestCubeName
    #   Move the MCG model cube to the destination
    os.system("mv "+MCGDict['CubeFile']+" "+FullDestName)
    #   Remove the MCG model
    os.system ("rm -r "+MCGDict['OutputName'])
        #   Finally adjust the header for the new file -- again, this is due to avoiding fish-eye effects
    AdjustHeaderPosBack(FullDestName,ObjDict['CubeFileName'])
    #   Similarly, there may be slight differences in the velocity reference point, so make sure they match -- note that the cubes do match regardless
    if SmoothSwitch == 0:
        AdjustHeaderVelBack(FullDestName,ObjDict['VelocityCubeFileName'])
    elif SmoothSwitch == 1:
        AdjustHeaderVelBack(FullDestName,ObjDict['SmoothedVelocityCubeFileName'])
    #   Add the provenance keywords to the cubes
    AddProvenanceKeyWords(FullDestName,FittingOptions['ProvenanceVals'],0,ObjDict)
    #   Add back the rotation keywords
    ReAddRotationKeys(FullDestName,ObjDict['CubeFileName'])
    #   Return the name of the MCG cubelet
    return FullDestName
    
def CopyFATFits(ObjDict,FolderDict,FittingOptions,AvgOutputDict,FATFncs):
    """
        This function copies the results from the individual FAT fits to the results directory so that an advanced user may explore them
    """
    #   Loop through all FAT fits
    for i in range(FittingOptions['nFATFits']):
        #   First make the specific FAT fits directory
        os.makedirs(AvgOutputDict['FATDir'][i], exist_ok=True)
        #   Get the name of the Parameter file
        FATFolder,ParamFile,CubeFile=FATFncs['GetFATNames'](ObjDict,FittingOptions['FATAnalysisFolders'][i],FolderDict)
        #   Adjust the name of the output file
        TargParamName=ParamFile.split(".")
        TargParamName=TargParamName[0].split("/")
        TargParamName=AvgOutputDict['OutputBaseName']+"_FATModel.txt"
    #   Now copy in the parameter and density files
        os.system("cp "+ ParamFile + " " + AvgOutputDict['FATDir'][i] + "/"+TargParamName)
        ConfigFile=FittingOptions['FATAnalysisFolders'][i]+FittingOptions['FATConfigFileName'][i]
        #   Do the same for the default configuration file name
        TargConfigName=FittingOptions['FATConfigFileName'][i].split(".")
        TargConfigName=AvgOutputDict['OutputBaseName']+"_FATInput.txt"
        #   And copy the config file to the FAT folder
        os.system("cp "+ ConfigFile + " " + AvgOutputDict['FATDir'][i] + "/"+TargConfigName)

def CopyBaroloFits(ObjDict,FolderDict,FittingOptions,AvgOutputDict,BaroloFncs):
    """
        This function copies the results for all the individual Barolo fits to the appropriate folders for expert users
    """
    #   Loop through all Barolo fits
    for i in range(FittingOptions['nBaroloFits']):
        #   First make the specific Barolo fits directory
        os.makedirs(AvgOutputDict['BaroloDir'][i], exist_ok=True)
        #   Get the name of the Parameter file
        BaroloFolder=BaroloFncs['NameBaroloAnalysisObjectFolder'](ObjDict,FittingOptions['BaroloAnalysisFolders'][i])
        ParamFile,DensFile,CubeFile=BaroloFncs['NameBaroloResultsFiles'](BaroloFolder,ObjDict)
        #   Adjust the target file names to work for ease of use
        NewParamName=AvgOutputDict['OutputBaseName']+"_BaroloModel.txt"
        NewDensName=AvgOutputDict['OutputBaseName']+"_BaroloSurfaceDensity.txt"
        #   Now copy in the parameter and density files
        os.system("cp "+ ParamFile + " " + AvgOutputDict['BaroloDir'][i] + "/"+NewParamName)
        os.system("cp "+ DensFile + " " + AvgOutputDict['BaroloDir'][i] + "/"+NewDensName)
        #   Finally copy in the input file
        InputFile=BaroloFncs['NameParamFile'](ObjDict,FolderDict,FittingOptions['BaroloFitNames'][i],0)
        #   As usual, adjust the name of the file to fit with the naming convention
        TargInFileName=AvgOutputDict['OutputBaseName']+"_BaroloInput.txt"
        InputFile=BaroloFolder+"/"+InputFile
        os.system("cp "+ InputFile + " " + AvgOutputDict['BaroloDir'][i] + "/"+TargInFileName)

def WriteAvgModel(AvgModel,ObjDict,AvgOutputDict,FittingOptions,FolderDict):
    """
        This function writes out a text file with kinematic model.
    """
    #   First get the date that the analysis is being run
    now = date.today()
    Date= str(now.year) + "-" + str(now.month) + "-" +  str(now.day)
    #   Add the date, source name, and model version to the dictionary for saving in the catalogue file
    AvgModel['Date']=Date
    AvgModel['SourceObs']=FolderDict['SourceName']
    AvgModel['ModelVers']=FittingOptions['ModelVersion']
    
    AvgModel['SBID']=FittingOptions['ProvenanceVals']['SBID']
    AvgModel['SRCTR']=FittingOptions['ProvenanceVals']['SRCTR']
    AvgModel['SRCVER']=FittingOptions['ProvenanceVals']['SRCVER']
    AvgModel['KINTR']=FittingOptions['ProvenanceVals']['KINTR']
    AvgModel['KINVER']=FittingOptions['ProvenanceVals']['KINVER']
    #   Name the output file
    FitName=AvgOutputDict['OutputDirName']+"/"+AvgOutputDict['OutputBaseName']+"_AvgModel.txt"
    #   Write out a the file header
    f=open(FitName,"w")
    HeaderStr="Object:\t "+ObjDict['ObjName']+"\n"
    HeaderStr+="Date: \t" + Date+"\n"
    #HeaderStr+="Version: \t" + str(FittingOptions['ModelVersion'])+"\n"
    HeaderStr+="SBID: \t" + str(AvgModel['SBID'])+"\n"
    HeaderStr+="SRCTR: \t" + str(AvgModel['SRCTR'])+"\n"
    HeaderStr+="SRCVER: \t" + str(AvgModel['SRCVER'])+"\n"
    HeaderStr+="KINTR: \t" + str(AvgModel['KINTR'])+"\n"
    HeaderStr+="KINVER: \t" + str(AvgModel['KINVER'])+"\n"
    HeaderStr+="Flag: \t"+str(AvgModel['FitFlags'])+"\n\n"
    #   Next write out the geometric parameters
    #       Start with a brief header
    HeaderStr+="Geometry Parameters\n"
    HeaderStr+="Param Name\t Value \t Error \n"
    f.write(HeaderStr)
    #   And now add each parameter to the file
    GeoStr="X_model (pixels) \t" +str(AvgModel['XCENTER'][0])+"\t\t"+str(AvgModel['XCENTER_ERR'][0])+"\n"
    GeoStr+="Y_model (pixels) \t" +str(AvgModel['YCENTER'][0])+"\t\t"+str(AvgModel['YCENTER_ERR'][0])+"\n"
    GeoStr+="RA_model (degrees) \t" +str(AvgModel['RA'][0])+"\t\t"+str(AvgModel['RA_ERR'])+"\n"
    GeoStr+="DEC_model (degrees) \t" +str(AvgModel['DEC'][0])+"\t\t"+str(AvgModel['DEC_ERR'])+"\n"
    GeoStr+="Inc_model (degrees) \t" +str(AvgModel['INCLINATION'][0])+"\t\t"+str(AvgModel['INCLINATION_ERR'][0])+"\n"
    GeoStr+="PA_model (degrees) \t" +str(AvgModel['POSITIONANGLE'][0])+"\t\t"+str(AvgModel['POSITIONANGLE_ERR'][0])+"\n"
    GeoStr+="PA_model,g (degrees) \t" +str(AvgModel['PA_GLOBAL'])+"\t\t"+str(AvgModel['PA_GLOBAL_ERR'])+"\n"
    GeoStr+="VSys_model (km/s) \t" +str(AvgModel['VSYS'][0])+"\t\t"+str(AvgModel['VSYS_ERR'][0])+"\n\n"
    f.write(GeoStr)
    #   Then do the rotation curve
    #       Add a brief header first
    HeaderStr="Rotation Curve\n"
    HeaderStr+="nR=\t"+str(len(AvgModel['R']))+"\n"
    HeaderStr+="Rad \t\t VROT_model \t e_VRot_model \t e_VROT_model,inc \n"
    HeaderStr+="('') \t\t (km/s) \t\t    (km/s) \t (km/s)"
    f.write(HeaderStr)
    #       Now start the profile string
    ProfileStr="\n"
    #       And loop through the radial profile adding the rotation curve
    for i in range(len(AvgModel['R'])):
        ProfileStr+=str(AvgModel['R'][i])+"\t\t"+str(AvgModel['VROT'][i])+"\t\t"+ str(AvgModel['VROT_ERR'][i])+"\t\t"+str(AvgModel['VROT_INC_ERR'][i])+"\n"
    f.write(ProfileStr)
    #   Finally do the same thing for the surface density
    #       Again write a small header
    HeaderStr="\nSurface Density Profile\n"
    HeaderStr+="nR=\t"+str(len(AvgModel['R_SD']))+"\n"
    HeaderStr+="Rad  \t\t SD_model \t\t e_SD_model \t SD_FO_model \t e_SD_FO_inc,model \n"
    HeaderStr+="('')\t\t (Msol/pc^2) \t (Msol/pc^2) \t (Msol/pc^2) \t (Msol/pc^2) \t "
    f.write(HeaderStr)
    #   Then go through the full surface density profile
    ProfileStr="\n"
    for i in range(len(AvgModel['R_SD'])):
        ProfileStr+=str(AvgModel['R_SD'][i])+"\t\t"+ str(AvgModel['SURFDENS'][i])+"\t\t"+str(AvgModel['SURFDENS_ERR'][i])+"\t\t"+str(AvgModel['SURFDENS_FACEON'][i])+"\t\t"+str(AvgModel['SURFDENS_INC_ERR'][i])+"\n"
    f.write(ProfileStr)
    #   Finally close the kinematic model file.
    f.close()
    
def WriteRCFile(AvgModel,ObjDict,AvgOutputDict,FittingOptions,FolderDict):
    """
        This function writes the rotation curve to a FITS file as a binary table
    """
    #   Name the output file
    FitName=AvgOutputDict['OutputDirName']+AvgOutputDict['OutputBaseName']+"_ModelRotationCurve.fits"
    #   Set up the FITS columns for the rotation curve
    col1 = fits.Column(name='Rad', format='E', array=AvgModel['R'],unit='arcsec')
    col2 = fits.Column(name='VRot_model', format='E', array=AvgModel['VROT'],unit='km/s')
    col3 = fits.Column(name='e_VRot_model', format='E', array=AvgModel['VROT_ERR'],unit='km/s')
    col4 = fits.Column(name='e_VRot_model,inc', format='E', array=AvgModel['VROT_INC_ERR'],unit='km/s')
    FitsTable=fits.BinTableHDU.from_columns([col1,col2,col3,col4])
    #   Make an empty primary hdu for the fits file
    hdu=fits.PrimaryHDU()
    #   Set some keywords in the primary header and table header for redundancy
    WriteBinTableHeaderKeyWords(hdu,AvgModel,ObjDict)
    WriteBinTableHeaderKeyWords(FitsTable,AvgModel,ObjDict)
    #   Combine the primary hdu and the table into an 'hdulist object'
    HduList=fits.HDUList([hdu,FitsTable])
    #   Write the Hdu file
    HduList.writeto(FitName,overwrite=True,output_verify='fix')
    #   Add the provenance keywords to the primary headers
    AddProvenanceKeyWords(FitName,FittingOptions['ProvenanceVals'],0,ObjDict)
    
    #   Add a fake image to the header for use with CASDA
    MakeFake1DHeader(FitName,ObjDict['CubeFileName'],ObjDict)
    
    
def WriteSDFile(AvgModel,ObjDict,AvgOutputDict,FittingOptions,FolderDict):
    """
        This function writes the surface density profile to a FITS file as a binary table
    """
    #   Name the output file
    FitName=AvgOutputDict['OutputDirName']+AvgOutputDict['OutputBaseName']+"_ModelSurfaceDensity.fits"
    #   Set up the surface density columns
    col1 = fits.Column(name='Rad_SD', format='E', array=AvgModel['R_SD'],unit='arcsec')
    col2 = fits.Column(name='SD_model', format='E', array=AvgModel['SURFDENS'],unit='M_sun/pc2')
    col3 = fits.Column(name='e_SD_model', format='E', array=AvgModel['SURFDENS_ERR'],unit='M_sun/pc2')
    col4 = fits.Column(name='SD_FO_model', format='E', array=AvgModel['SURFDENS_FACEON'],unit='M_sun/pc2')
    col5 = fits.Column(name='e_SD_FO_model,inc', format='E', array=AvgModel['SURFDENS_FACEON_ERR'],unit='M_sun/pc2')
    #   Write the columns to a FITS binary table
    FitsTable=fits.BinTableHDU.from_columns([col1,col2,col3,col4,col5])
    #   Make an empty primary hdu for the fits file
    hdu=fits.PrimaryHDU()
    #   Set some keywords in the primary header and table header for redundancy
    WriteBinTableHeaderKeyWords(hdu,AvgModel,ObjDict)
    WriteBinTableHeaderKeyWords(FitsTable,AvgModel,ObjDict)
    #   Combine the primary hdu and the table into an 'hdulist object'
    HduList=fits.HDUList([hdu,FitsTable])
    #   Write the Hdu file
    HduList.writeto(FitName,overwrite=True,output_verify='fix')
    
    #   Add the provenance keywords to the primary headers
    AddProvenanceKeyWords(FitName,FittingOptions['ProvenanceVals'],0,ObjDict)
    
    #   Add a fake image to the header for use with CASDA
    MakeFake1DHeader(FitName,ObjDict['CubeFileName'],ObjDict)
    
def WriteGeoFile(AvgModel,ObjDict,AvgOutputDict,FittingOptions,FolderDict):
    """
        This function writes the geometric parameters to a FITS file as a binary table
    """
    #   Name the output file
    FitName=AvgOutputDict['OutputDirName']+AvgOutputDict['OutputBaseName']+"_ModelGeometry.fits"
    #   Make an array of the geometric parameter names
    NameArr=['X_model','Y_model','RA_model','DEC_model','Inc_model','PA_model','PA_model,g','VSys_model']
    #   Make an array of the geometric parameter keywords used in the Tilted Ring dictionary
    KeyArr=['XCENTER','YCENTER','RA','DEC','INCLINATION','POSITIONANGLE','PA_GLOBAL','VSYS']
    #   Make an array of the units for the various parameters
    UnitsArr=np.array(['pixels','pixels','degrees','degrees','degrees','degrees','degrees','km/s'])
    #   Set up an empty array for the measurements, uncertainties, and names
    ValArr=np.array([None]*len(NameArr))
    ErrArr=np.array([None]*len(NameArr))
    TestNameArr=np.array([None]*len(NameArr))
    #   Start a counter for the array
    i=0
    #   Loop through all the keywords
    for x in NameArr:
        #   Set the current keyword, associated error, and parameter name
        key=KeyArr[i]
        keyerr=key+"_ERR"
        TestNameArr[i]=x
        #   Check if it's an array or not
        if isinstance(AvgModel[key],np.ndarray) or isinstance(AvgModel[key],list):
            #   If it is, store the first value (as a flat disk, it should be constant
            ValArr[i]=AvgModel[key][0]
        else:
            #   If not, store the value itself
            ValArr[i]=AvgModel[key]
        #   Do the same for the uncertainties
        if isinstance(AvgModel[keyerr],np.ndarray):
            ErrArr[i]=AvgModel[keyerr][0]
        else:
            ErrArr[i]=AvgModel[keyerr]
        #   Update the counter
        i+=1
    #   Put all the arrays into columns
    col1 = fits.Column(name='ParamName', format='A10', array=TestNameArr)
    col2 = fits.Column(name='Units', format='A8', array=UnitsArr)
    col3 = fits.Column(name='Value', format='E', array=ValArr)
    col4 = fits.Column(name='Error', format='E', array=ErrArr)
    #   Make the fits table
    FitsTable=fits.BinTableHDU.from_columns([col1,col2,col3,col4])
    #   Make an empty primary hdu for the fits file
    hdu=fits.PrimaryHDU()
    #   Set some keywords in the primary header and table header for redundancy
    WriteBinTableHeaderKeyWords(hdu,AvgModel,ObjDict)
    WriteBinTableHeaderKeyWords(FitsTable,AvgModel,ObjDict)
    #   Combine the primary hdu and the table into an 'hdulist object'
    HduList=fits.HDUList([hdu,FitsTable])

    #   Write the Hdu file
    HduList.writeto(FitName,overwrite=True,output_verify='fix')
    
    #   Add the provenance keywords to the primary headers
    AddProvenanceKeyWords(FitName,FittingOptions['ProvenanceVals'],0,ObjDict)
    #   Add a fake image to the header for use with CASDA
    MakeFake1DHeader(FitName,ObjDict['CubeFileName'],ObjDict)

    
def WriteBinTableHeaderKeyWords(hdu,AvgModel,ObjDict):

    hdu.header.set('OBJECT',ObjDict['ObjName'])
    hdu.header.set('DATE',AvgModel['Date'])
    hdu.header.set('ORIGIN','WKAPP')

def WriteOutputTable(AvgModels,FittingOptions):
    """
        This function saves the table of model parameters to a single csv file
    """
    #   Start by setting up the header
    HeaderStr="Preliminary WALLABY Kinematic Models\n"
    #   Next make the PANDAS dataframe
    DF=MakePANDASDataFrame(AvgModels)
    #   Name the CSV file
    csvname=FittingOptions['AvgModelReleaseBaseName']+".csv"
    #   Write the dataframe to the csv file
    DF.to_csv(csvname,index=False)
    #   Place the csv file into the appropriate folder
    MoveScript="mv "+csvname+" "+FittingOptions['AverageModelFolder']+"/."
    os.system(MoveScript)

def MakePANDASDataFrame(AvgModels):
    """
        This function converts the list of TR dictionaries into a single PANDAS dataframe
    """
    #   Start with an empty dictionary
    FullDict={}
    #   Figure out the number of models
    nModels=np.shape(AvgModels)[0]
    #   Make an empty list for every parameter to be stored in the dataframe
    IncList=[None]*nModels
    IncErrList=[None]*nModels
    PAList=[None]*nModels
    PAErrList=[None]*nModels
    PA_GList=[None]*nModels
    PA_GErrList=[None]*nModels
    VSysList=[None]*nModels
    VSysErrList=[None]*nModels
    XList=[None]*nModels
    XErrList=[None]*nModels
    YList=[None]*nModels
    YErrList=[None]*nModels
    NameList=[None]*nModels
    FolderList=[None]*nModels
    IDList=[None]*nModels
    FlagsList=[None]*nModels
    DateList=[None]*nModels
    SourceObsList=[None]*nModels
    VersList=[None]*nModels
    RAList=[None]*nModels
    RAErrList=[None]*nModels
    DECList=[None]*nModels
    DECErrList=[None]*nModels
    SBIDList=[None]*nModels
    SRCTRList=[None]*nModels
    SRCVERList=[None]*nModels
    KINTRList=[None]*nModels
    KINVERList=[None]*nModels
    ra_sofiaList=[None]*nModels
    dec_sofiaList=[None]*nModels
    freq_sofiaList=[None]*nModels
    #   Loop through all the models in the field
    for i in range(nModels):
        #   Add the model values to the appropriate lists
        IncList[i]=AvgModels[i]['INCLINATION'][0]
        IncErrList[i]=AvgModels[i]['INCLINATION_ERR'][0]
        PAList[i]=AvgModels[i]['POSITIONANGLE'][0]
        PAErrList[i]=AvgModels[i]['POSITIONANGLE_ERR'][0]
        PA_GList[i]=AvgModels[i]['PA_GLOBAL']
        PA_GErrList[i]=AvgModels[i]['PA_GLOBAL_ERR']
        VSysList[i]=AvgModels[i]['VSYS'][0]
        VSysErrList[i]=AvgModels[i]['VSYS_ERR'][0]
        XList[i]=AvgModels[i]['XCENTER'][0]
        XErrList[i]=AvgModels[i]['XCENTER_ERR'][0]
        YList[i]=AvgModels[i]['YCENTER'][0]
        YErrList[i]=AvgModels[i]['YCENTER_ERR'][0]
        RAList[i]=AvgModels[i]['RA'][0]
        RAErrList[i]=AvgModels[i]['RA_ERR']
        DECList[i]=AvgModels[i]['DEC'][0]
        DECErrList[i]=AvgModels[i]['DEC_ERR']
        IDList[i]=AvgModels[i]['ID']
        NameList[i]=AvgModels[i]['NAME']
        FolderList[i]=AvgModels[i]['FOLDER']
        DateList[i]=AvgModels[i]['Date']
        SourceObsList[i]=AvgModels[i]['SourceObs']
        VersList[i]=AvgModels[i]['ModelVers']
        FlagsList[i]=AvgModels[i]['FitFlags']
        SBIDList[i]=AvgModels[i]['SBID']
        SRCTRList[i]=AvgModels[i]['SRCTR']
        SRCVERList[i]=AvgModels[i]['SRCVER']
        KINTRList[i]=AvgModels[i]['KINTR']
        KINVERList[i]=AvgModels[i]['KINVER']
        ra_sofiaList[i]=AvgModels[i]['ra_sofia']
        dec_sofiaList[i]=AvgModels[i]['dec_sofia']
        freq_sofiaList[i]=AvgModels[i]['freq_sofia']
        
    #   Save all the lists to the full dictionary
    FullDict={'name':NameList,'ra':ra_sofiaList,'dec':dec_sofiaList,'freq':freq_sofiaList,'team_release':SRCTRList,'team_release_kin':KINTRList,'Vsys_model':VSysList,'e_Vsys_model':VSysErrList, 'X_model':XList,'e_X_model':XErrList,'Y_model':YList,'e_Y_model':YErrList,'RA_model':RAList,'e_RA_model':RAErrList,'DEC_model':DECList,'e_DEC_model':DECErrList,'Inc_model':IncList,'e_Inc_model':IncErrList,'PA_model':PAList,'e_PA_model':PAErrList,'PA_model,g':PA_GList,'e_PA_model,g':PA_GErrList,'QFlag_model':FlagsList}
    
    #       Now set up the csv to store the profile.
    #   First set the profile keys for saving to the header
    ProfileKeys=['Rad','Vrot_model','e_Vrot_model','e_Vrot_model,inc','Rad_SD','SD_model','e_SD_model','SD_FO_model','e_SD_FO_model,inc']
    #   Next set the corresponding keys from the average model dictionary
    ModelKeys=['R','VROT','VROT_ERR','VROT_INC_ERR','R_SD','SURFDENS','SURFDENS_ERR','SURFDENS_FACEON','SURFDENS_INC_ERR']
    #   Initialize the arrays in the full dictionary
    for x in ProfileKeys:
        FullDict[x]=np.array([None]*nModels)
    #   Loop through all the models
    for i in range(nModels):
        j=0
        #   For each model, loop through all the keys
        for x in ProfileKeys:
            #   Get the model key
            y=ModelKeys[j]
            #   Select the array
            Arr=AvgModels[i][y]
            #   Turn the array into a string with commas to work with CASDA
            ArrStr = ','.join([str(elem) for elem in Arr])
            #   Store the string
            FullDict[x][i]=ArrStr
            #   Increment the counter
            j+=1
    
    #   Convert the dictionary to a dataframe
    DF=pd.DataFrame.from_dict(FullDict)
    return DF

def AdjustHeaderPosBack(CubeName,OriCubeName):
    """
        This function adjusts a cube header back to the original SoFiA header.  This is done to deal with 'fish-eye' effects.
    """
    #   Open the current cube
    Cube=fits.open(CubeName,mode='update')
    #   Open the original cube
    OriCube=fits.open(OriCubeName)
    #   Set the list of spatial keywords that need adjusting
    Keys=['CRPIX1','CRVAL1','CRPIX2','CRVAL2']
    #   Set the current cube keywords to the original cube's keywords
    for x in Keys:
        Cube[0].header.set(x,OriCube[0].header[x])
    #   Also adjust the angular unit to deg instead of DEGREES for consistency
    Cube[0].header.set('CUNIT1','deg')
    Cube[0].header.set('CUNIT2','deg')
    #   Close the original cube
    OriCube.close()
    #   Save the current cube
    Cube.flush()
    #   Close the current cube
    Cube.close()
    
def AdjustHeaderVelBack(CubeName,OriCubeName):
    """
        This function matches the velocity keywords in a cube header to some original cube.
    """
    #   Open the current cube
    Cube=fits.open(CubeName,mode='update')
    #   Open the original cube
    OriCube=fits.open(OriCubeName)
    #   Set the list of velocity keywords that need adjusting
    Keys=['CRPIX3','CRVAL3']
    #   Set the current cube keywords to the original cube's keywords
    for x in Keys:
        Cube[0].header.set(x,OriCube[0].header[x])
    #   Close the original cube
    OriCube.close()
    #   Save the current cube
    Cube.flush()
    #   Close the current cube
    Cube.close()
    
def ReAddRotationKeys(CubeName,OriCubeName):
    """
        This function matches the velocity keywords in a cube header to some original cube.
    """
    #   Open the current cube
    Cube=fits.open(CubeName,mode='update')
    #   Open the original cube
    OriCube=fits.open(OriCubeName)
    #   Set the list of velocity keywords that need adjusting
    Keys=['LONPOLE','LATPOLE']
    #   Set the current cube keywords to the original cube's keywords
    for x in Keys:
        Cube[0].header.set(x,OriCube[0].header[x])
    #   Close the original cube
    OriCube.close()
    #   Save the current cube
    Cube.flush()
    #   Close the current cube
    Cube.close()


def AddProvenanceKeyWords(CubeName,ProvDict,HeaderID,ObjDict):
    """
        This function adds the provenance key words to a header using the provenance dictionary
    """
    #   Open the current cube
    Cube=fits.open(CubeName,mode='update')
    #   Set the current cube keywords to the original cube's keywords
    for x in ProvDict.keys():
        Cube[HeaderID].header.set(x,ProvDict[x])
    #   Also adjust the object key word to the correct name
    Cube[HeaderID].header.set('OBJECT',ObjDict['ObjName'])
    #   And make sure the date is correct
    now = date.today()
    Date= str(now.year) + "-" + str(now.month) + "-" +  str(now.day)
    Cube[HeaderID].header.set('DATE',Date)
    #   Save the current cube
    Cube.flush()
    #   Close the current cube
    Cube.close()

def MakeFake1DHeader(CubeName,OriCubeName,ObjDict):
    #   Open the current cube
    Cube=fits.open(CubeName,mode='update')
    #   Add a 1 pixel array for a 'fake' image needed for CASDA
    Cube[0].data = np.array([[[0]]]).astype("int16")
    #   Open the original cube
    OriCube=fits.open(OriCubeName)
    WCSLIBAxes=3
    CRPIX1=0
    CRPIX2=0
    CRPIX3=0
    CUNIT1="deg"
    CUNIT2="deg"
    CUNIT3="Hz"
    #OBJECT=ObjDict['ObjName']
    CRVAL1=ObjDict['CatEntry']['ra']
    CRVAL2=ObjDict['CatEntry']['dec']
    CRVAL3=ObjDict['CatEntry']['freq']
  
    PredefinedDict={'WCSAXES':WCSLIBAxes,'CRPIX1':CRPIX1,'CRPIX2':CRPIX2,'CRPIX3':CRPIX3,'CUNIT1':CUNIT1,'CUNIT2':CUNIT2,'CUNIT3':CUNIT3,'CRVAL1':CRVAL1,'CRVAL2':CRVAL2,'CRVAL3':CRVAL3}
    for x in PredefinedDict.keys():
        Cube[0].header.set(x,PredefinedDict[x])
        
    OriKeys=['CDELT1','CDELT2','CDELT3','CTYPE1','CTYPE2','CTYPE3']
    for x in OriKeys:
        Cube[0].header.set(x,OriCube[0].header[x])
    #   Close the original cube
    OriCube.close()
    #   Save the current cube
    Cube.flush()
    #   Close the current cube
    Cube.close()
