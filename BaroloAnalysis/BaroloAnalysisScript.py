import numpy as np
import astropy
import copy as copy
import os

from . import BaroloConfigOptions as BCO
from . import BaroloModelAnalysis as BMA

"""
    This module contains routines for running the Barolo analysis on the WALLABY  galaxies.  It contains the routines:
    BaroloAnalysisFolderPrep --> Set up a folder for the Barolo Analysis results
    RunBaroloScript --> Run Barolo for some galaxy
    WriteBaroloParamFile --> Write the Barolo parameter file that controls the Barolo run
    SetBaroloPreambleStr --> Write the preamble to the Barolo parameter file
    WriteFittingOptions --> Write the fitting options to the Barolo parameter file
    AddGeometricAvgValues --> Add the averaged geometric values to the Barolo parameter file
    AverageGeometry --> Average the geometric values from a Barolo run
    AdjustPA --> Adjust the position angle to be between 0-360 degrees
    AvgArr --> Average an array
    AvgVSys -->  Get the average for the Vsys from the 1st run, with accounting for bad rings.
    LoadFirstRunFit --> Do a bsic load in the results from a Barolo run
"""


def BaroloAnalysisFolderPrep(BaroloOptions):
    """
        This function makes a folder for each Barolo run to be done.
    """
    #   Loop through all Barolo fitting options and set them
    for i in range(BaroloOptions['nFits']):
        os.makedirs(BaroloOptions['BaroloAnalysisFolders'][i], exist_ok=True)
        



def RunBaroloScript(ObjDict,FolderDict,BaroloOptions):
    """
        This function performs a flat disk estimate using 3D-Barolo for a specific galaxy.
        Due to issues with the Barolo 2-stage mode for making flat disks, we've implemented our own '2-stage' method that fixes the geometric parameters in the second run.
    """
    #   First grab the Barolo fitting options set in BaroloConfigOptions.py
    FitSwitches=BCO.FittingOptions()
    
    #   It may be necessary to run Barolo on different sets of cubes of different resolutions, so loop over all the total number of possible fits.
    for i in range(BaroloOptions['nFits']):
        #   Make an empty array for the first barolo fits for passing purposes
        FirstBaroloFit=[]
        #   Get the name of the output folder
        OutputFolder=BMA.NameBaroloAnalysisObjectFolder(ObjDict,BaroloOptions['BaroloAnalysisFolders'][i])
        #   Run Barolo twice for each object
        for j in range(2):
            #   Write the Barolo Parameter File
            ParamFileName=WriteBaroloParamFile(ObjDict,FolderDict,BaroloOptions,OutputFolder,FitSwitches,i,j,FirstBaroloFit)
            #   Run Barolo using the written parameter file
            BaroloCmd=FolderDict['BaroloPath']+ " -p " + ParamFileName
            os.system(BaroloCmd)
            #   Move the parameter file into the target directory
            BaroloCleanCmd="mv "+ParamFileName + " " + OutputFolder +"/."
            os.system(BaroloCleanCmd)
            #   If this is the fit step, load in the Barolo Model
            if j ==0:
                FirstBaroloFit,FailedFit=LoadFirstRunFit(OutputFolder)
                #       If the first fit has failed, don't run again
                if FailedFit:
                    print("First fit failed")
                    break
def WriteBaroloParamFile(ObjDict,FolderDict,BaroloOptions,OutputFolder,FitSwitches,FitStep,RunStep,FirstFitParams):
        """
            This function writes the parameter file that controls the 3D-Barolo analysis.
        """
        #   Get the preamble for the Barolo parameter file
        BaroloScriptStr=SetBaroloPreambleStr(ObjDict,BaroloOptions,OutputFolder)
        #   Add the Mask options to the Barolo scipt
        print(BaroloOptions)
        BaroloScriptStr+=BCO.SetBaroloMaskingOptions(ObjDict,BaroloOptions['UseSoFiAMasks'])
        #   Set the fitting options
        BaroloScriptStr+=WriteFittingOptions(FitSwitches,RunStep)
        #   If this is the first run and there are prior estimates, these need to added to the script
        if RunStep==0:
            if BaroloOptions['FitNames'][FitStep] == 'Optical':
                BaroloScriptStr+=AddPriorEstimates(ObjDict)
        #   If this is the second run, then we need to add the geometric values to the file
        if RunStep==1:
            BaroloScriptStr+=AddGeometricAvgValues(FirstFitParams)
    
        #   Once the string is fully put together, write out the parameter file
            #   First get the parameter file name
        ParamFileName=BMA.NameParamFile(ObjDict,FolderDict,BaroloOptions['FitNames'][FitStep],RunStep)
        ParamFile=open(ParamFileName,'w')
        ParamFile.write(BaroloScriptStr)
        ParamFile.close()
        return ParamFileName
    
def SetBaroloPreambleStr(ObjDict,BaroloOptions,OutputFolder):
    """
        Write the initial portion of the Barolo parameter file that will be needed regardless of the specific fitting options
    """
    #   Figure out if using the normal cubes or the smoothed data
    if BaroloOptions['UseSmoothedCubes']:
        TargCube=ObjDict['SmoothedVelocityCubeFileName']
    else:
        TargCube=ObjDict['VelocityCubeFileName']
        
    #   Add the name of the FITS cube to be analyzed
    PreambleStr="FITSFILE\t"+TargCube+"\n"
    #   Write out the folder to place the analysis into
    PreambleStr+="OUTFOLDER\t"+OutputFolder+"\n"
    #   Set the 3D fitting options to true
    PreambleStr+="3DFIT\t true\n"
    #   Set the linear parameter
    PreambleStr+="LINEAR 0.42\n"
    #   Only use a single thread as these are small cubelets.
    PreambleStr+="\nTHREADS 1\n"
    return PreambleStr
    
def WriteFittingOptions(FitSwitches,RunStep):
    """
        This function writes out the specific fitting options as specified in the BaroloConfigOptions.py file
    """
    #   Initialize an empty string
    FitStr="\n"
    #   Loop through all options in the Fitting options (FitSwitches) dictionary
    for x in FitSwitches:
        #   For the parameters that are fit (i.e. Free to vary), the analysis will depend on which step in the 2-step fitting that we've implemented
        if x == 'FREE':
            #   For the 2nd run, only fit the rotation curve.
            if RunStep ==1:
                FitStr+="FREE\tVROT\n"
                #   Otherwise use the parameters listed in BaroloConfigOptions.py
            else:
                FitStr+=x+"\t"+FitSwitches[x]+"\n"
        #   For all other parameters, write down the name of the switch and the string associated with it.
        else:
            FitStr+=x+"\t"+FitSwitches[x]+"\n"
    return FitStr

    
def AddGeometricAvgValues(FirstFitParams):
    """
        On the 2nd Barolo run, the geometric parameters should be fixed to the average values found in the first run.
    """
    #   Initially average all the geometric parameters in the first run
    AvgFit=AverageGeometry(FirstFitParams)
    
    #   Then write out a geometry string that fixes the values to the average values.  These are used as initial estimates in Barolo.  Since those parameters are not free to vary (down in the WriteFittingOptions function here), they remain constant in the 2nd run.
    GeoValsStr="\n"
    GeoValsStr+="VSYS\t"+str(AvgFit['VSys_Avg']) +"\n"
    GeoValsStr+="XPOS\t"+str(AvgFit['X_Avg']) +"\n"
    GeoValsStr+="YPOS\t"+str(AvgFit['Y_Avg']) +"\n"
    GeoValsStr+="INC\t"+str(AvgFit['Inc_Avg']) +"\n"
    GeoValsStr+="PA\t"+str(AvgFit['PA_Avg']) +"\n"
    return GeoValsStr
    
def AverageGeometry(FirstFitParams):
    """
        This function averages the geometric parameters from Barolo fit.  It takes in the array of 'first fit params' that were loaded in from the model file.
    """
    print("Averaging Geometry")
    #   Get the geometric parameters
    Inc=FirstFitParams[4]
    PA=FirstFitParams[5]
    X=FirstFitParams[9]
    Y=FirstFitParams[10]
    VSys=FirstFitParams[11]
    #   Adjust the PA to be between 0-360, and deal with PA's near 0
    PA=AdjustPA(PA)
    #   Average the parameters
    FitDict={}
    #   VSys needs special treatment due to potential bad rings
    FitDict['VSys_Avg']=AvgVSys(VSys)
    FitDict['Inc_Avg']=AvgArr(Inc)
    FitDict['X_Avg']=AvgArr(X)
    FitDict['Y_Avg']=AvgArr(Y)
    FitDict['PA_Avg']=AvgArr(PA)
    #   If the average PA is above or below 360 due to the cyclic adjustment, fix the average
    FitDict['PA_Avg']=FitDict['PA_Avg']%360
    return FitDict

def AdjustPA(PA):
    """
        Occasionally, Barolo has problems with a the position angle of a ring or 2 where it is off by some number of arbitrary rotations.  This breaks the inbuilt 2-stage averaging.  This function fixes those position angles so that our averaging will not be affected by such problems.
    """
    #   Adjust the Position angle for averaging purposes
    #   If there are more than a single ring in the Barolo fit, adjust each PA
    try:
        #   Loop through each ring
        for i in range(len(PA)):
            #   First get the angle between 0-360 degrees
            PA[i]=PA[i]%360.
        #   If PA is around 0 degrees, there might be an issue in the averaging process due to the cyclical nature of the PA
        #   To deal with this, loop through the position angles
        for i in range(len(PA)):
            #   Check the difference between this ring and the center ring
            Diff=PA[i]-PA[0]
            #   If the difference is more than 180 degrees, the averaging will be affected, so we need to adjust the angles such that they are less than 180 degrees -- this typically occurs when the PA ~0
            if np.abs(Diff) >=180.:
                #   For negative differences, adjust the ring's PA upwards.
                if Diff <0.:
                    PA[i]=PA[i]+360.
                #   For positive differences, adjust the ring's PA downwards.
                else:
                    PA[i]=PA[i]-360.
    #   If it's only a single ring, there are no arrays to deal with.
    except:
        PA=PA%360
    return PA
    
def AvgArr(Arr):
    """
        This function averages an array or returns the specific value if it's a float.
    """
    #   Try doing the array operations
    try:
        #   Sum up all elements of the array
        Tot=np.sum(Arr)
        #   Average them by dividing the sum by the length of the array.
        Avg=Tot/float(len(Arr))
    #   If only a float, return that
    except:
        Avg=Arr
    return Avg
    
def AvgVSys(Arr):
    """
        The Vsys average requires special treatment as it has the potential to go much farther from the 'true' value in the outermost rings.
    """
    #   First average the array using all values
    VSysAvg_Ini=AvgArr(Arr)
    #   Check that the array is more than one element, as this issue can only happen in such a case
    print(type(Arr))
    print(Arr)
    if isinstance(Arr, (np.float64)):
        return VSysAvg_Ini
    #   Set the acceptance threshold
    Threshold=20.
    #   Make an empty list to start
    VSysTemp=[]
    #   Get the temporary array and it's size
    VSysTemp,count=MakeThresholdArr(Arr,VSysAvg_Ini,Threshold)
    #  Make sure the temporary array is large enough for the averaging
    while count<=1:
        Threshold=Threshold+20
        VSysTemp,count=MakeThresholdArr(Arr,VSysAvg_Ini,Threshold)
        #print("Threshold Problem",count,Threshold,VSysTemp)
    #   Average the new array
    VSysAvg=AvgArr(VSysTemp)
    Avg=VSysAvg
    return Avg
    
def MakeThresholdArr(Arr,Avg,Threshold):
    """
        This function makes a new array of values within some threshold of the average
    """
    #   Set the new array and a counter
    NewArr=[]
    count=0
    #   Loop through the array
    for i in range(len(Arr)):
        #   Get the difference from the average
        Diff=Arr[i]-Avg
        #   Check that the difference is below the acceptable threshold
        if np.abs(Diff) < Threshold:
            #   If it is, append the value to the temporary array
            NewArr.append(Arr[i])
            count+=1
    #   Return the new array and count
    return NewArr,count
    
def LoadFirstRunFit(OutputFolder):
    """
        This function loads in the results from the 1st Barolo analysis as a basic numpy array
    """
    #   Get the name of the results file -- Barolo always names this #path/rings_final1.txt
    ResultsFile=OutputFolder+"/rings_final1.txt"
    #   Try to open the results of the fit if it was successful
    try:
        #   Load in the fits
        FirstRunResults=np.loadtxt(ResultsFile,skiprows=1,unpack=True)
        #   Change the name of the file so that it's kept for later runs
        os.system("mv " + ResultsFile + " " + OutputFolder+"/ring_final1_1stpass.txt")
        FailedFit=False
    #   If the file doesn't load, set the FailedFit flag to exit the Barolo analysis.
        if np.shape(FirstRunResults)[0] ==0:
            FailedFit=True
    except:
        FailedFit=True
        FirstRunResults=[]
    return FirstRunResults,FailedFit
