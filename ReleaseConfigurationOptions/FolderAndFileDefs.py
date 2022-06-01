from . import HydraField as HyaF
from . import HydraDR1_Field as HyaFDR1
from . import NormaField as NorF
from . import NGC4636Field as N4636F
"""
    This file contains routines that define file/folder locations for the kinematic analysis.  The routines are:
    BasicFolderAndFileLoc --> Sets the basic file/folder definitions
    ExcelCatDetails --> Sets the details of the SoFiA excel sheet
    InitializeObjectDictionary --> Sets properties/definitions for a specific galaxy
"""


def BasicFolderAndFileLoc():
    """
        Set the basic file/folder definitions for a given field
    """
    
    ParentFolder="/Users/nate/Dropbox/WALLABY/"
    SizeLimit=2.
    
    #       The Hydra field options
    FolderDict=HyaF.FolderAndFileLoc(ParentFolder)
    #       The Norma field options
    #FolderDict=NorF.FolderAndFileLoc(ParentFolder)
    #       The NGC4636 field options
    #FolderDict=N4636F.FolderAndFileLoc(ParentFolder)
    #       The Hydra field options
    #FolderDict=HyaFDR1.FolderAndFileLoc(ParentFolder)
    
    FolderDict['ParentFolder']=ParentFolder
    FolderDict['SizeLimit']=SizeLimit
    
    return FolderDict


def ExcelCatDetails():
    """
        Set the details for loading in the SoFiA catalogue as an excel file for a specific field
    """
    #   The Hydra field excel sheet options
    ExcelDetailsDict=HyaF.ExcelCatDetails()
    #   The Norma field excel sheet options
    #ExcelDetailsDict=NorF.ExcelCatDetails()
        #   The NGC4636 field excel sheet options
    #ExcelDetailsDict=N4636F.ExcelCatDetails()
    #   The Hydra field excel sheet options
    #ExcelDetailsDict=HyaFDR1.ExcelCatDetails()
    
    return ExcelDetailsDict
    
    

def InitializeObjectDictionary(WallabyCat,ID,FolderDict):
    """
        This function sets up up an object dictionary that contains variables/definitions specific to a given galaxy.  The galaxy itself is selected by the ID and uses information in both the SoFiA WALLABY catalogue and the folder dictionary to set the values.
    """
    #   First get the index from the SoFiA Catalogue by matching the ID
    RowIndx=WallabyCat.index[WallabyCat['id']==ID][0]
    #   Select the entire row of SoFiA values
    WallabyCatValues=WallabyCat.loc[RowIndx]
    
    #   Get the name of the object from the catalogue name
    Name=WallabyCat['name'][RowIndx]
    #   Get the string Jx-y from the Name by splitting at the space
    NumStr=Name.split(' ')[1]
    #   Set up the name with an underscore as this is used in most naming conventions instead of spaces
    ObjName_Underscored="WALLABY_"+NumStr
    #   Get the base name/location for SoFiA files
    ObjBaseFileName=FolderDict['SofiADataFolder']+ObjName_Underscored
    #   Set both the SoFiA cubelet and mask names
    CubeName=ObjBaseFileName+"_cube.fits.gz"
    MaskName=ObjBaseFileName+"_mask.fits.gz"
    """
    #       The Hydra DR1 files use a different naming convention, so adjust them here
    ObjBaseFileName=FolderDict['SofiADataFolder']+FolderDict['BaseName']+"_"+NumStr
    CubeName=ObjBaseFileName+"_cube.fits"
    MaskName=ObjBaseFileName+"_mask.fits"
    """
    #   Set the names of both the smoothed and unsmoothed velocity cubelets.  The specific files are generated in preprocessing sets (WallabyPreliminary.py)
    VelocityCubeName=FolderDict['VelFolder']+FolderDict['BaseName']+"_"+NumStr+"_cube_Vel.fits"
    SmoothedVelocityCubeName=FolderDict['SmoothedVelFolder']+FolderDict['BaseName']+"_"+NumStr+"_cube_Vel_h3.fits"
    
    
    #   Initialize a dicitonary for any parameter estimates from non-kinematic modelling sources
    NonKinematicParameterEstimates={}
    #   Place all the information into a dictionary
    GalaxyDict={'RowIndx':RowIndx,'CatEntry':WallabyCatValues,'ObjName':Name,'NumStr':NumStr,'ObjName_Underscored':ObjName_Underscored,'ObjFileBaseName':ObjBaseFileName,'CubeFileName':CubeName,'MaskFileName':MaskName,'VelocityCubeFileName':VelocityCubeName,'SmoothedVelocityCubeFileName':SmoothedVelocityCubeName,'NonKinematicEstimates':NonKinematicParameterEstimates}
    return GalaxyDict
