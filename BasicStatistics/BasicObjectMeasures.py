import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from spectral_cube import SpectralCube
from astropy import wcs

import copy as copy
import os

from . import AstroFncs as AF

"""
    This module contains routines for making most of the basic measurements of the galaxies in a particular SoFiA catalogue.  It contains the routines:
    SizeEstimate --> Estimates the size of galaxies.
    DistanceEstimate --> Estimates the distance to galaxies based on their redshift
    MassCalc --> Estimates the HI mass in a detection
    RHI_Estimate --> Estimates the size of the galaxy based on the estimated mass
    VelElements --> Estimates the number of velocity elements
    CubeVolEstimate--> Estimates the total channel volume of a cube
    SNEstimate --> Measures the S/N in a cubelet in a variety of ways
"""


def SizeEstimate(WallabyCat):
    """
        This function gives the size of galaxies in the SoFiA catalogue using the Ell_maj parameter.  It currently has the pixelsize and beamsize hardcoded into the results.
    """
    #   Set the pixel and beam size as well as the pixel/beam ratio
    pixelSize=1.66666666667E-03
    beamSize=8.33333333333E-03
    PixelBeamRatio=beamSize/pixelSize
    #   The SoFiA Ell_maj column is in pixels, so convert them to beams by the pixel/beam ratio.
    ObjectSize_Beams=WallabyCat['ell_maj']/PixelBeamRatio
    return ObjectSize_Beams
    
def DistanceEstimate(WallabyCat):
    """
        This function estimates the distance to a galaxy based on the redshift.  It relies on a hardcoded value for H0 and a rest frequency.  It is worth noting that this value is not used in any other portions of the proto-pipeline and is included only for testing purposes.  The distance estimates are not reported as part of WKAPP.
    """
    #   Set H0 and the 21 cm rest frequency
    H0=70
    RestFreq=1.42040575179E+09
    #   Select the central frequency of all SoFiA detections from the catalogue.
    CatFreq=WallabyCat['freq']
    #   Get the redshift velocity from these frequencies in km/s.
    CatVel=AF.RedShiftConv(WallabyCat['freq'],RestFreq)/1000.
    #   Get the distance using the Hubble flow.
    Distance=AF.DistEst(H0,CatVel)
    return Distance

def MassCalc(WallabyCat):
    """
        Estimate the total HI mass in each detection.  This depends on the distance estimate (and therefore also H0).  As with the distances, this is not used further in the proto-pipeline, nor is it reported as part of WKAPP.  It has mostly been used for testing purposes during code development.
    """
    #   Set the frequency bin size.
    dF=  WallabyCat['df']
    #   Get the total mass using the total flux, distance estimate, and bin size.
    Mass=AF.MassFromFlux(WallabyCat['Distance'],WallabyCat['f_sum'],dF)
    return Mass

def RHI_Estimate(WallabyCat):
    """
        Estimate the value of R_HI from the HI mass-size relation.  As witht he distance and mass estimates (on which this estimate relies), it is done for testing purposes only and is not included in the final WKAPP analysis.
    """
    #   Calculate the value of R_HI in kpc based on the estimated HI mass.
    RHI=AF.CalcRHI(WallabyCat['logM'])
    #   Convert the size of R_HI from kpc to arcsec-->Note that this should cancel out the dependence on the distance estimate.
    RHI_Projected=AF.ProjectedSize_ArcSecond(RHI,WallabyCat['Distance'])
    return RHI, RHI_Projected

def VelElements(WallabyCat):
    """
        This estimates the number of velocity channels covered by a detection using the catalogue w20 measurement and bin size.  This was used during code development for some tests, but is not used anywere else.
    """
    nVelElements=WallabyCat['w20']/WallabyCat['df']
    return nVelElements

def CubeVolEstimate(WallabyCat):
    """
        This estimates the number of independant spaxil elements covered by a detection using the size of the object in beams and velocity channels.  This was used during code development for some tests, but is not used anywere else.  It gives a volume in beams^2*channels.
    """
    VolEstimate=WallabyCat['ObjectSize_Beams']**2.*WallabyCat['VelElements']
    return VolEstimate




def SNEstimate(ObjDict,WallabyCat):
    """
        This function estimates the S/N for a particular galaxy in a number of different ways.  It returns the peak S/N, average S/N, and integrated S/N (based on Eq. A1 in Westmeier et al. 2021).
    """
    #   Start with the cube rms from the SoFiA catalogue
    RMS=WallabyCat['rms'][ObjDict['RowIndx']]
    #   To get the signal, it is necessary to mask the data cube.
    #       Start by setting the mask cube name
    MaskName=ObjDict['MaskFileName']
    #   Open up the data cube.
    DC=fits.open(ObjDict['CubeFileName'])
    #   Open up the mask cube file
    Mask=fits.open(MaskName)
    #   Mask the data cube
    MaskedData=DC[0].data*Mask[0].data
    
    #   Get the peak flux inside the masked data cube.
    MaxSig=np.nanmax(MaskedData)
    #   Get the peak S/N by dividing the
    SN_Peak=MaxSig/RMS
    #   Calculate the average signal in each cell that is inside the mask
    AvgSig=np.nansum(MaskedData)/np.nansum(Mask[0].data)
    #   Get the average S/N by dividing the average signal by the RMS per cell from the SoFiA analysis.
    SN_Avg=AvgSig/RMS
    #   Get the integrated (or observed) S/N based on Eq. A1 in Westmeier et al. 2021.
    #       This measurement includes the area of the beam in pixels in the denominator, so calculate that first.
    Beam_Pixels=DC[0].header['BMAJ']/np.abs(DC[0].header['CDELT1'])
    BeamArea=(2.*np.pi/2.355**2.)*Beam_Pixels**2.
    #   Now apply Eq. A1 by summing up the masked signal and dividing it by the square root of the number of cells multiplied by the beam area in pixels.  Finally divide the whole thing by the RMS per cell to get a unitless S/N.
    SN_Obs=np.nansum(MaskedData)/np.sqrt(np.nansum(Mask[0].data)*BeamArea)/RMS

    #   Close the data cube and the mask cube.
    DC.close()
    Mask.close()
    #   Return the RMS, and different S/N measures for use in other portions of WKAPP.  In particular SN_Obs is often used as a selecting a subsample for kinematic analysis.
    return RMS,SN_Peak,SN_Avg,SN_Obs
