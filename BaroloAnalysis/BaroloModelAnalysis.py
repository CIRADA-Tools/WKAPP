import numpy as np
import astropy
import copy as copy
import os

"""
    This module contains routines for analysing the results from a 3D-Barolo analysis.  The main use for these routines is loading in the results from a Barolo analysis into a tilted ring dictionary.  This dictionary is then used to make diagnostic plots and generate an average model with other fits.  It contains the functions:
    BaroloModelFnc -->  This function stores the other functions into a dictionary for use in other modules without directly importing this module.
    NameParamFile --> This function gets the name of the Barolo input parameter file.  This is mostly used when organizing and saving the final models.
    NameBaroloAnalysisObjectFolder --> This function gets the name of the folder that holds the specific galaxy's Barolo analysis results.
    NameBaroloResultsFiles --> This function gets the name of the Barolo results files (parameters and surface density).
    LoadBaroloModel --> This function loads in the Barolo results from an analysis and stores it into a TR dictionary.
    LoadBaroloParamFile --> This function loads in most of the Barolo model parameters.
    LoadBaroloSurfDens --> This function loads in the Barolo surface density profile.
    BaroloParamstoDict --> This function takes the Barolo parameters and converts them to a dictionary.
    BadBaroloFit --> This cuntion is used when Barolo has failed for some reason to make an empty dictionary and flag the result as an issue.
"""


def BaroloModelFnc():
    """
        This function stores the other routines into a dictionary for easy use without importing the netire module.
    """
    BaroloFncDict={'NameBaroloAnalysisObjectFolder':NameBaroloAnalysisObjectFolder,'LoadBaroloModel':LoadBaroloModel,'NameBaroloResultsFiles':NameBaroloResultsFiles,'NameParamFile':NameParamFile}
    return BaroloFncDict

def NameParamFile(ObjDict,FolderDict,FitName,RunStep):
    """
        This function gets the name of the specific parameter file.  This is mostly used for saving files once the average model is generated.
    """
    #   First get the base name
    PartialName=FolderDict['BaseName']+"_"+ObjDict['NumStr']+"_"+FitName
    #   Change the name of the file depending on which run step this is
    if RunStep==0:
        ParamFile=PartialName+".param"
    else:
        ParamFile=PartialName+"_2ndPass.param"
    return ParamFile

def NameBaroloAnalysisObjectFolder(ObjDict,SpecificBaroloFolder):
    """
        This module gets the folder name that holds the results of the Barolo analysis for a specific galaxy.
    """
    FolderName=SpecificBaroloFolder+ObjDict['NumStr']
    return FolderName
    
def NameBaroloResultsFiles(ResultsFolder,ObjDict):
    """
        This function gets the names of the output model files from the Barolo analysis.
    """
    #   Get the name of the file with the Barolo model
    ParamFile=ResultsFolder+"/rings_final1.txt"
    #   Get the name of the file with the Barolo surface density
    DensFile=ResultsFolder+"/densprof.txt"
    #   Get the name of the Barolo model cube.
    CubeFile=ResultsFolder+"/"+ObjDict['ObjName_From_Cube']+"mod_azim.fits"
    return ParamFile,DensFile,CubeFile
    

def LoadBaroloModel(ObjDict,FolderDict,FitFolder,FitLabel):
    """
        Load in all the results of a Barolo analysis into a dictionary.
    """
    #   First get the folder name for the Barolo model
    ResultsFolder=NameBaroloAnalysisObjectFolder(ObjDict,FitFolder)
    #   Get the names of the specific Barolo results files
    ResultsFile,DensFile,CubeFile=NameBaroloResultsFiles(ResultsFolder,ObjDict)
    print("Barolo results filename", ResultsFile, os.path.isfile(ResultsFile))
    #   Load in the parameter file
    BaroloParams,FitAchieved=LoadBaroloParamFile(ResultsFile)
    #   Next try to load in the surface density profile...if the model file was loaded in correctly
    if FitAchieved:
        SDParams,FitAchieved=LoadBaroloSurfDens(DensFile)
    #   If the fit failed, set the bad fit dictionary
    if FitAchieved ==False:
        BaroloDict=BadBaroloFit()
    #   If the fit was successful, then there is a number of steps that need to happen
    else:
        #   First turn the model parameter arrays into a dictionary
        BaroloDict=BaroloParamstoDict(BaroloParams,SDParams)
        #   Turn the switch that this fit can be used in the model averaging algorithm
        BaroloDict['FitForAveraging']=True
    #   Set the Barolo Dict to have the label for the particular fit
    BaroloDict['Label']=FitLabel
    return BaroloDict

def LoadBaroloParamFile(ParamFile):
    """
        Load in the results from a Barolo output into a numpy array
    """
    #   Start by assuming the fit is true.
    FitAchieved=True
    #   Try to load in the model file.
    try:
        Params=np.loadtxt(ParamFile,skiprows=1)
    except:
        #   If the fit fails, set the fit flag to false and make an empty parameter array
        FitAchieved=False
        Params=[]
    #   Occasionally Barolo makes a file with only the header but no rings.  Check for that by looking at the shape of the parameters array.  If it's empty also set the succesful fit flag to false.
    if np.shape(Params)[0] == 0:
        FitAchieved=False
    return Params,FitAchieved

def LoadBaroloSurfDens(DensFile):
    """
        Load in the surface density profile from a Barolo run into a numpy array.
    """
    #   Try to load in the surface density profile
    try:
        SDParams=np.loadtxt(DensFile,skiprows=13,usecols=(0,9,10,11),unpack='True')
        FitAchieved=True
        #   If the load fails, set the successful fit flag to false and make an empty surface density array.
    except:
        FitAchieved=False
        SDParams=[]
    return SDParams,FitAchieved



def BaroloParamstoDict(Params,SDParams):
    """
        This routine converts the numpy model and surface density arrays loaded in from a Barolo analysis into a consistent tilted ring dictionary that can be used with other tilted ring modelling ring results.
    """

    #   First set the core dictionary from the parameter arrays.  The specific columns are selected from the Barolo results file.
    TRParamsDict={'R':Params[:,1],'R_SD':SDParams[0,:],'XCENTER':Params[:,9] \
            ,'YCENTER':Params[:,10], 'INCLINATION':Params[:,4] \
            ,'POSITIONANGLE':Params[:,5], 'VSYS':Params[:,11], 'VROT':Params[:,2]\
            ,'VRAD':Params[:,12], 'VDISPERSION':Params[:,3],'Z0':Params[:,7]\
                ,'SURFDENS':SDParams[1,:],'SURFDENS_ERR':SDParams[2,:],'FITAchieved':True}
    #   The face-on SD is needed for generating MCG models later on
    #       However the face on SD must be converted to MCG units (M_sol/pc^2).  This needs the conversion factor between these 2 units (see )
    SDConv=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    TRParamsDict['SURFDENS_FACEON']=SDParams[3,:]*SDConv
    TRParamsDict['SURFDENS_FACEON_ERR']=SDParams[2,:]
    
    #   For some diagnostics, Barolo is run with uncertainties on, but this is not always true.  So try to add the error terms to the dictionary if they exist.  If they are not present, make an array of zeros as some functions assume that the arrays exist.
        #       Get a rotation error if it exists.
    try:
        VerrLow=Params[:,13]
        VerrHigh=Params[:,14]
    except:
        VerrLow=Params[:,1]*0.
        VerrHigh=Params[:,1]*0.
        #       Get the inclination error if it exists
    try:
        IncerrLow=Params[:,15]
        IncerrHigh=Params[:,16]
    except:
        IncerrLow=Params[:,1]*0.
        IncerrHigh=Params[:,1]*0.
        #       Get the position angle error if it exists
    try:
        PAerrLow=Params[:,17]
        PAerrHigh=Params[:,18]
    except:
        PAerrLow=Params[:,1]*0.
        PAerrHigh=Params[:,1]*0.

    #       Place the errors into the TR dictionary
    VErr=[VerrLow,VerrHigh]
    TRParamsDict['VROT_ERR']=VErr
    IncErr=[IncerrLow,IncerrHigh]
    TRParamsDict['INCLINATION_ERR']=IncErr
    PAErr=[PAerrLow,PAerrHigh]
    TRParamsDict['POSITIONANGLE_ERR']=PAErr
    return TRParamsDict

def BadBaroloFit():
    """
        Even when Barolo fails, a model dictionary needs to exist in other functions.  So here, make an empty dictionary with flags set to false as necessary.
    """
    Disk={'R':[],'XCENTER':[] \
    ,'YCENTER':[], 'INCLINATION':[] \
    ,'POSITIONANGLE':[], 'VSYS':[], 'VROT':[]\
    ,'VRAD':[], 'VDISPERSION':[],'Z0':[]\
    ,'SURFDENS':[],'FITAchieved':False,'RESID':[]
    ,'FitForAveraging':False,'CHI2':[]}
    return Disk
