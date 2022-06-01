import numpy as np


import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

import copy as copy
import os

"""
    This file contains routines for converting SoFiA cubelets to a more usable form for FAT and 3DBarolo.  It contains the routines:
    SetUpVelFolders --> Sets up folders to hold both the 4 km/s and 12 km/s velocity cubelets.
    CubeHeaderConvert --> Converts a SoFiA cubelet in frequency space to velocity space and recenters the spatial reference point for the cubelets.  It saves these cubelets to a the full-resolution velocity cubelet folder.
    SetNewSpatialRefCoordinates --> Sets the spatial CRPIX and CRVALs to the center of the cube.
    Frequency_VelocityConversion --> Sets the spectral header values to velocity units.
"""

def SetUpVelFolders(FolderDict):
    """
        This function makes a folder to store the 4 km/s velocity cubelets ('VelFolder') and 12 km/s velocity cubelets ('SmoothedVelFolder').  The names of these folders must be specified in the FolderDict dictionary.
    """
    #   Check if the target velocity folder exists -- if not, make it
    os.makedirs(FolderDict['VelFolder'], exist_ok=True)
    #   Do the same for the smoothed velocity cube folder
    os.makedirs(FolderDict['SmoothedVelFolder'], exist_ok=True)
    
    
def CubeHeaderConvert(ObjDict,FolderDict,AstroFnc):
    """
        This function adjusts the header of the SoFiA cubelets so that the output cubes are in velocity space with a reference point in the center of the cubelet.
    """
    #   Get the target name for the converted velocity file.
    OutName=ObjDict['VelocityCubeFileName']
    #   Open up the frequency file
    Cubehdu=fits.open(ObjDict['CubeFileName'])
    #   Read the frequency cube header
    CubeHeader=Cubehdu[0].header
    #   Get the wcs coordinate information for the frequency cube
    w = wcs.WCS(Cubehdu[0].header)
    #   First do the spatial conversion
    Cubehdu=SetNewSpatialRefCoordinates(Cubehdu,w,CubeHeader,AstroFnc)
    #   Now to the frequency-velocity conversion
    Cubehdu=Frequency_VelocityConversion(Cubehdu,CubeHeader,AstroFnc)
    #   Write out the cube with the new header to the given file name
    try:
        Cubehdu[0].writeto(OutName)
    except:
        os.system("rm " + OutName)
        Cubehdu[0].writeto(OutName)
    #   Close the file
    Cubehdu.close()


def SetNewSpatialRefCoordinates(CubeHDU,WCSObj,CubeHeader,AstroFnc):
    """
        This function sets the spatial reference point of some cubelet to the center of the cubelet.  It is built on the wcslib function pix2world
    """
    #   First specify the new reference pixels to be the middle of the cube
    RefChannel=int(CubeHeader['NAXIS3']/2)  #   This is needed for the pix2world function
    RefX=int(CubeHeader['NAXIS1']/2)
    RefY=int(CubeHeader['NAXIS2']/2)
    
    #   Set the new reference pixel as a proper array
    NewReferencePix=[[RefX,RefY,RefChannel]]
    #   Get the coordinates for this pixel using the wcslib
    NewReferenceCoords=WCSObj.wcs_pix2world(NewReferencePix,1)
    
    #   Save the new reference position to the cube header
    CubeHDU[0].header.set('CRPIX1',NewReferencePix[0][0])
    CubeHDU[0].header.set('CRVAL1',NewReferenceCoords[0,0])
    CubeHDU[0].header.set('CRPIX2',NewReferencePix[0][1])
    CubeHDU[0].header.set('CRVAL2',NewReferenceCoords[0,1])

    #   Also set the cube units to degrees --> this is needed for many other codes
    CubeHDU[0].header.set('CUNIT1','DEGREE')
    CubeHDU[0].header.set('CUNIT2','DEGREE')
    #   return the cube header
    return CubeHDU

def Frequency_VelocityConversion(CubeHDU,CubeHeader,AstroFnc):
    """
        This function converts the spectral axis from frequency to velocity.  It assumes that the units of the frequency axis are Hertz.
        
        Because the source cubes are so large, the size of the velocity channels must be calculated at the central reference point of the cube and not at the overall reference point as dV is not constant across the larger cube.
    """
    #      First get the frequecy cube reference channel, frequency, and bin size
    RefChannel=CubeHDU[0].header['CRPIX3']
    RefFrequency=CubeHDU[0].header['CRVAL3']
    DeltaFrequency=CubeHDU[0].header['CDELT3']
    
    #   Now set the new reference channel to be the middle of the cube
    NewRefChannel=int(CubeHeader['NAXIS3']/2)

    #   Get the number of channels being changed
    DeltaChannel=NewRefChannel-RefChannel

    #   Get the frequency at the new reference channel
    F0=RefFrequency+(DeltaChannel)*DeltaFrequency
    #   Get the frequency at the channel right after the new reference channel
    F1=RefFrequency+(DeltaChannel+1)*DeltaFrequency
    
    #   Set the rest frequency to that of the HI 21 cm line in Hz
    RestFreq=1.42040575179E+09
    #   Convert the reference frequency to a velocity in m/s
    V0=AstroFnc['RedShiftConv'](F0,RestFreq)
    #   Convert the frequency in the neighbouring channgel to a velocity in m/s
    V1=AstroFnc['RedShiftConv'](F1,RestFreq)
    #   Use these two velocities to get the velocity width for channels in the cubelet
    dV=V1-V0

    #   Save all the new velocity information to the CubeHDU
    CubeHDU[0].header.set('CRPIX3',NewRefChannel)
    CubeHDU[0].header.set('CRVAL3',V0)
    CubeHDU[0].header.set('CDELT3',dV)
    CubeHDU[0].header.set('CTYPE3','VOPT')
    CubeHDU[0].header.set('CUNIT3','m/s')
    #   return the Cube HDU
    return CubeHDU
