import numpy as np
import astropy
from astropy.io import fits
import copy as copy
import os

"""
    This file contains routines to run the 3DBarolo smoothing routines on the velocity converted cubelets to make ones at 12 km/s resolution.  It does this by writing a 3DBarolo input file and lets Barolo handle the smoothing process.  This file contains the routines:
    BaroloSmoothing --> Smooths the cubelets with a Hanning window of 3 channels using the 3DBarolo smoothing functionality.
"""

def BaroloSmoothing(ObjDict,FileDict):
    """
        This function smooths cubelets to 3x lower resolution.  For WALLABY it goes from 4 km/s to 12 km/s resolution.  It writes the appropriate Barolo input file and then runs the code.
    """
    #   Set the path to the 3DBarolo executable
    BaroloPath=FileDict['BaroloPath']

    #Write the barolo parameter file necessary to smooth the cubes
    TempScript=WriteBaroloSmoothingInputFile(ObjDict,FileDict)
    #   Run the 3DBarolo in the smoothing mode
    BaroloScript=BaroloPath+" -p "+ TempScript
    os.system(BaroloScript)
    #   Remove the temporary script
    os.system("rm "+TempScript)
    #   Change the Barolo output name to a slightly more appropriate name
    NewName=ObjDict['SmoothedVelocityCubeFileName']
    #   Adjust the header to clean up the object keyword
    VelHDU=fits.open(ObjDict['VelocityCubeFileName'])
    try:
        HDUObjName=VelHDU[0].header['OBJECT']
        OutNameBase=HDUObjName.replace(' ','_')
    except:
        OutNameBase='NONE'
    VelHDU.close()
    #   Move the Barolo file to the smoothed function with the appropriate naming convention
    BaroloOutputName=FileDict['SmoothedVelFolder']+OutNameBase+"_h3_red.fits"
    MoveScript="mv "+BaroloOutputName+ " "+NewName
    os.system(MoveScript)

def WriteBaroloSmoothingInputFile(ObjDict,FileDict):
    #Write the barolo parameter file necessary to smooth the cubes
    BaroloFile="FITSFILE\t"+ObjDict['VelocityCubeFileName']+"\n"
    BaroloFile+="SMOOTHSPEC \t true\n"
    BaroloFile+="WINDOW_TYPE\ t HANNING\n"
    BaroloFile+="WINDOW_SIZE \t 3\n"
    BaroloFile+="REDUCE \t true\n"
    BaroloFile+="OUTFOLDER \t"+FileDict['SmoothedVelFolder'] +"\n"

    #   Write the script to a temporary file
    TempScript="TempBaroloSmoothScript.par"
    g=open(TempScript,"w")
    for x in BaroloFile:
        g.write(x)
    g.close()
    
    return TempScript
