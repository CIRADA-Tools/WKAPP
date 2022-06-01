import numpy as np
import astropy
import copy as copy
import os

"""
    This module contains routines for loading in the results of a FAT analysis for a galaxy and storing them into a pair of tilted ring dictionaries.  The reason for using two dictionaries is that, even in a flat disk fit, FAT separates the galaxy into an approaching and receding side.  It contains the routines:
    FATModelFnc --> This routine stores the FAT routines that might be called elsewhere as a dictionary.
    LoadFATModel --> This function loads in the results of FAT analysis and stores it into a pair of dictionaries.
    GetFATNames --> This routine gets the names of the FAT results file and model cube
    ReadFATParamsFile --> This routine reads the specifc FAT results file.
    FATErrAdjust --> This routine places the uncertainties from the two disks into common arrays.
    ExtractLineValues --> This routine gets the values in a specific line of the FAT results file.
    ReadFATKeyWord --> This routine reads the FAT keyword in a given line.
    BadFATFit -->  This routine is used when the FAT analysis fails.  It makes an empty dictionary and sets the successful fit flag to false.
    
"""

def FATModelFnc():
    """
        This routine places useful FAT routines into a dictionary so that other modules can use them without explicitly importing this module.
    """
    FATFncDict={'LoadFATModel':LoadFATModel,'GetFATNames':GetFATNames}
    return FATFncDict

def LoadFATModel(ObjDict,FolderDict,FitFolder,FitLabel):
    """
        This module loads in the results from a FAT analysis of a specific galaxy and places them into a pair of dictionaries.
    """
    #       First get the names of the FAT results folder, file, and cube file
    FATFolder,ParamsFile,CubeFile=GetFATNames(ObjDict,FitFolder,FolderDict)
    print("FAT Params file", ParamsFile)
    #   Next read in the FAT parameters file
    Disk1,Disk2=ReadFATParamsFile(ParamsFile)
    
    #   If the FAT run was successful set a few extra dictionary elements for both comparison with other runs and averaging together into the a combined model.
    if Disk1['FITAchieved']:
        #   FAT only gives a face-on surface density, so the SD dictionary elements need to be set properly.
        #       Since velocity and SD are on the same radial grid (different than 3DBarolo), set a surface density radial grid
        Disk1['R_SD']=Disk1['R']
        #   Adjust the errors for the Disk1 dictionary.
        Disk1=FATErrAdjust(Disk1)
        #   Flag the fit for use in the averaging process.
        Disk1['FitForAveraging']=True
        Disk1['SURFDENS_FACEON']=Disk1['SURFDENS']  #The surface density is already face-on, but needs to be set for MCG later on
    #   Do the same set of things for the second dictionary
    if Disk2['FITAchieved']:
        Disk2['R_SD']=Disk2['R']
        Disk2=FATErrAdjust(Disk2)
        #   Since the second disk is a copy of the first (in flat disk models), do not use it when averaging different fits together.
        Disk2['FitForAveraging']=False
        Disk2['SURFDENS_FACEON']=Disk2['SURFDENS']  #The surface density is already face-on, but needs to be set for MCG later on
    #   Set the disk1 and disk2 labels from those set in the field configuration options (see routines in ReleaseConfigurationOptions/)
    Disk1['Label']=FitLabel[0]
    Disk2['Label']=FitLabel[1]
    return Disk1,Disk2
    
def GetFATNames(ObjDict,FitFolder,FolderDict):
    """
        This function gets the names of the specific galaxy analysis folder, the results file, and the model cube.
    """
    FATFolder=FitFolder+FolderDict['BaseName']+"_"+ObjDict['NumStr']
    ParamsFile=FATFolder+"/Finalmodel/FinalModel.def"
    CubeFile=FATFolder+"Finalmodel/FinalModel.fits"
    return FATFolder,ParamsFile,CubeFile

def ReadFATParamsFile(FileName):
    """
        This routine reads a FAT results file and places the results into a pair of dictionaries for each disk.
    """
    #   Try to open the final file
    try:
        f=open(FileName,"r")
        #   If the file can't be opened, note that there is no model and set up the failed fit dictionaries.
    except:
        print("No Final FAT Model")
        Disk1=BadFATFit()
        Disk2=BadFATFit()
        return Disk1,Disk2
    #   Break the results file into individual lines.
    FullFile=f.readlines()  #   Read file into lines
    #   Set up a set of characters that are used to separate values in a specific line
    SpaceChar=(" ","   ","  ")
    #   Since the file exists, initialize the dictionaries and set the fit achieved flag to true.
    Disk1={'FITAchieved':True}
    Disk2={'FITAchieved':True}
    #   In a FAT flat disk file, the model parameters are listed in lines 17-55
    for i in range(17,55,1):
        #   For each line, get the specific value/variable and update the appropriate dictionary
        Disk1,Disk2=ExtractLineValues(FullFile[i],SpaceChar,Disk1,Disk2)
    #   MCG does need a set of radial velocities, so set them to zero in the dictionaries.
    Disk1['VRAD']=Disk1['R']*0.
    Disk2['VRAD']=Disk1['R']*0.
    #   The 2nd disk also needs a radial grid in the dictionary.
    Disk2['R']=Disk1['R']
    #   Close the results file.
    f.close()
    return Disk1,Disk2
    
    
def FATErrAdjust(Disk):
    """
        This routine adjusts the uncertainties in the FAT TR dictionary to ones that better fit the overall sructure.  The odd structures are needed in the diagnostic plots.
    """
    Disk['VROT_ERR']=[Disk['VROT_ERR'],Disk['VROT_ERR']]
    Disk['POSITIONANGLE_ERR']=[Disk['PA_ERR'],Disk['PA_ERR']]
    Disk['INCLINATION_ERR']=[Disk['INC_ERR'],Disk['INC_ERR']]
    Disk['SURFDENS_ERR']=[Disk['R']*0.,Disk['R']*0.]
    return Disk

def ExtractLineValues(LineStr,Spacer,FATDict1,FATDict2):
    """
        This routine takes a single line from a FAT results file and figures out the keyword/parameter that the line refers to and all the measurements for that parameter.
    """
    #       First split the line by the equals sign
    EqualSplit=LineStr.split("=")
    #       Get the keyword and the line formating
    ReadFormat,VarName,DiskSwitch=ReadFATKeyWord(EqualSplit[0])
    if ReadFormat == -1:
        #       If it's not one of the recognized keywords, don't do anything
        return FATDict1,FATDict2
    else:
        #       Remove the \n character from the line
        values=EqualSplit[1].split("\n")[0]
        #       Split into an array according to the spacer and format for the variable
        Arr=values.split()
    #   Turn the array into floats
    ArrFloat=np.array(list(map(float,Arr)))
    #   Add the array to the Disk1 or Disk2 dictionaries by the keyword
    if DiskSwitch==0:
        FATDict1[VarName]=ArrFloat
    else:
        FATDict2[VarName]=ArrFloat
    return FATDict1,FATDict2

def ReadFATKeyWord(KeyStr):
    """
        This routine figures out the dictionary keyword from a FAT key string and whether it is for the Disk! or Disk2 dictionaries.
    """
    #   First get rid of the leading # sign
    TestStr=KeyStr.strip('#')
    #   Next get rid of any remaining whitespaces
    TestStr=TestStr.strip()
    #   Set an empty string, disk switch and read format in case no valid words are found
    VarName=[]
    DiskSwitch=0
    #   Check all the Disk1 keywords
    ReadFormat=-1
    if TestStr == 'NUR':
        ReadFormat=0
        VarName='nR'
    elif TestStr == 'RADI':
        ReadFormat=0
        VarName='R'
    elif TestStr == 'VROT':
        ReadFormat=0
        VarName='VROT'
    elif TestStr == 'VROT_ERR':
        ReadFormat=2
        VarName='VROT_ERR'
    elif TestStr == 'Z0':
        ReadFormat=0
        VarName='Z0'
    elif TestStr == 'SBR':
        ReadFormat=0
        VarName='SURFDENS'
    elif TestStr == 'INCL':
        ReadFormat=0
        VarName='INCLINATION'
    elif TestStr == 'INCL_ERR':
        ReadFormat=2
        VarName='INC_ERR'
    elif TestStr == 'PA':
        ReadFormat=0
        VarName='POSITIONANGLE'
    elif TestStr == 'PA_ERR':
        ReadFormat=2
        VarName='PA_ERR'
    elif TestStr == 'XPOS':
        ReadFormat=0
        VarName='RA'
    elif TestStr == 'YPOS':
        ReadFormat=0
        VarName='DEC'
    elif TestStr == 'VSYS':
        ReadFormat=0
        VarName='VSYS'
    elif TestStr == 'SDIS':
        ReadFormat=0
        VarName='VDISPERSION'
    elif TestStr == 'SDIS_ERR':
        ReadFormat=2
        VarName='VDISP_ERR'
    #   Check all the Disk2 keywords.
    elif TestStr == 'VROT_2':
        ReadFormat=0
        VarName='VROT'
        DiskSwitch=1
    elif TestStr == 'VROT_2_ERR':
        ReadFormat=2
        VarName='VROT_ERR'
        DiskSwitch=1
    elif TestStr == 'Z0_2':
        ReadFormat=0
        VarName='Z0'
        DiskSwitch=1
    elif TestStr == 'SBR_2':
        ReadFormat=0
        VarName='SURFDENS'
        DiskSwitch=1
    elif TestStr == 'INCL_2':
        ReadFormat=0
        VarName='INCLINATION'
        DiskSwitch=1
    elif TestStr == 'INCL_2_ERR':
        ReadFormat=2
        VarName='INC_ERR'
        DiskSwitch=1
    elif TestStr == 'PA_2':
        ReadFormat=0
        VarName='POSITIONANGLE'
        DiskSwitch=1
    elif TestStr == 'PA_2_ERR':
        ReadFormat=2
        VarName='PA_ERR'
        DiskSwitch=1
    elif TestStr == 'XPOS_2':
        ReadFormat=0
        VarName='RA'
        DiskSwitch=1
    elif TestStr == 'YPOS_2':
        ReadFormat=0
        VarName='DEC'
        DiskSwitch=1
    elif TestStr == 'VSYS_2':
        ReadFormat=0
        VarName='VSYS'
        DiskSwitch=1
    elif TestStr == 'SDIS_2':
        ReadFormat=0
        VarName='VDISPERSION'
        DiskSwitch=1
    elif TestStr == 'SDIS_2_ERR':
        ReadFormat=2
        VarName='VDISP_ERR'
        DiskSwitch=1
    return ReadFormat,VarName,DiskSwitch

def BadFATFit():
    """
        Even when FAT fails, a model dictionary needs to exist in other functions.  So here, make an empty dictionary with flags set to false as necessary.
    """
    Disk={'R':[],'XCENTER':[] \
    ,'YCENTER':[], 'INCLINATION':[] \
    ,'POSITIONANGLE':[], 'VSYS':[], 'VROT':[]\
    ,'VRAD':[], 'VDISPERSION':[],'Z0':[]\
    ,'SURFDENS':[],'FITAchieved':False,'FitForAveraging':False}
    return Disk
