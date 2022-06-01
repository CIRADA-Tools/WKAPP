import numpy as np
import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

"""
    This module contains a number of useful and simple astronomy functions that are used throughout WKAPP.  The specific functions in this module are:
    FuncDict --> This places all the other functions in this module into a dictionary for use in other routines without directly importing this module.
    RedShiftConv --> This converts a redshift to a velocity.
    DistEst --> This estimates a distance based on the Hubble flow.
    MassFomFlux --> This calculates a total HI mass from a total flux in Jy
    CalcRHI --> This estimates the size of a galaxy from the HI mass.
    ProjectedSize_ArcSecond --> This calculates a size in arcseconds from kpc.
    CalcRA_Dec_FromCube_And_Center --> This gets the central RA and DEC from a pixel point
    CalcCenter_FromCube_And_RADEC --> This gets a pixel point from a pair of RA and DEC.
"""


def FuncDict():
    """
        Place useful functions in the rest of this module into a dictionary.
    """
    AstronFncDict={'RedShiftConv':RedShiftConv,'DistEst':DistEst,'MassFromFlux':MassFromFlux,'CalcRHI':CalcRHI,'ProjectedSize_ArcSecond':ProjectedSize_ArcSecond,'CalcRA_Dec_FromCube':CalcRA_Dec_FromCube_And_Center,'CalcCenter_FromCube':CalcCenter_FromCube_And_RADEC}
    return AstronFncDict

def RedShiftConv(Freq,RestFreq):
    """
        Calculate the velocity of some object by it's redshifted spectrum
    """
    #   First get the redshift
    z=(RestFreq-Freq)/Freq
    #   Next get the velocity by multiplying the redshift by the speed of light in m/s.
    Vel=z*2.9979245e8
    return Vel
    
def DistEst(H0,Vel):
    """
        Get the distance of some galaxy using the Hubble flow by dividing the velocity in km/s by a supplied H0.
    """
    Distance=Vel/H0
    return Distance
    
def MassFromFlux(Distance,Flux,dF):
    """
        Get the HI mass from the flux using the () relation.  It assumes that the flux is in units of Jy Hz
    """
    #   This assumes that the flux is in units of Jy freq
    Mass=0.236*(Distance*1000.)**2.*Flux/abs(dF)
    #   Store the mass in logarithmic units
    Mass=np.log10(Mass.astype(float))
    return Mass
   
def CalcRHI(logMHI):
    """
        Use the HI mass-size relation to estimate R_HI
    """
    #Use the HI mass - diameter relation from
    #Wang+2016 to estimate RHI from logMHI.
    #-->INPUT: log(MHI/Mo)
    #-->OUTPUT: RHI in kpc
    #--------------------
    #   Set the slope and intercept of the mass-size relation
    slope = 0.506
    intercept = -3.293
    #   Get the logarithmic diameter
    logDHI = slope*logMHI + intercept
    #   Convert to a linear radius
    RHI = (10**logDHI)/2.
    return RHI
   
   
def ProjectedSize_ArcSecond(Size,Distance):
    """
        This function simply converts a size in kpc to arcseconds. Note that this assumes the size is in kpc and the distance is in Mpc.
    """
    Size_P=Size/(Distance*1000) * 206265
    return Size_P



def CalcRA_Dec_FromCube_And_Center(X,Y,CubeDict):
    """
        This function uses the astropy WCSLIB functionality to convert a point in pixels to a value in RA-DEC coordinates.  It requires the cube header for the conversion.
    """
    #   Set the wcslib functionality and cube header infor to the variable w
    w=CubeDict['CubeWCS']
    #   Make an empty list of RA and Dec
    RA=[]
    DEC=[]
    #   Loop through all pairs of X and Y
    for i in range(np.shape(X)[0]):
        #   Set a 3D pixel coordinate, but we are not interested in the spectral coordinate.
        pixcrd=[[X[i], Y[i],0]]
        #   Use the wcslib functions to convert the pixel coordinates to equatorial coordinates.
        CentCoord=w.wcs_pix2world(pixcrd,0)
        #   Append the RA and DEC values to the existing list.
        RA.append(CentCoord[0,0])
        DEC.append(CentCoord[0,1])

    return RA,DEC

def CalcCenter_FromCube_And_RADEC(RA,DEC,CubeDict):
    """
        This function uses the astropy WCSLIB functionality to convert a point in RA-DEC coordinates to pixels.  It requires the cube header for the conversion.
    """
    #   Set the wcslib functionality and cube header infor to the variable w
    w=CubeDict['CubeWCS']
    #   Make an empty list of pixel coordinates
    X=[]
    Y=[]
    #   Loop through all pairs of RA and DEC
    for i in range(len(RA)):
            #   Set a 3D real space coordinate, but we are not interested in the spectral coordinate.
        RealCoords=[[RA[i],DEC[i],0]]
        #   Use the wcslib functions to convert the equatorial coordinates to pixel coordinates.
        PixCoords=w.wcs_world2pix(RealCoords,0)
        #   Append the X and Y values to the existing list.
        X.append(PixCoords[0,0])
        Y.append(PixCoords[0,1])

    return X,Y
