import numpy as np
import astropy

"""
    This module contains routines that specify the options for running 3D-Barolo.  It contains the routines:
        FittingOptions --> This gives the specific fitting options for 3D-Barolo
        SetBaroloMaskingOptions --> This sets how the mask that Barolo uses will be generated.
"""


def FittingOptions():
    """
        This function sets the specific fitting options for 3D-Barolo and places them into a dictionary
    """
    #   Set the type of normalization to use
    Normalization="AZIM\n"
    #   Set the default velocity dispersion in km/s.
    DefaultDispersion="10"
    #   Set the Error estimate flag -- since this is not used in the averaging, false is a good choice.
    Errors="false"
    #   Set whether Barolo will do a 2 stage fit -- this is set false due to some small issues in Barolo's 2-stage fitting.
    TwoStage="false"
    #   Decide whether to make the default Barolo plots -- in later steps we use our own plotting scripts.
    PlotSwitch="false"
    #   Set the default radial seperation in arcsec -- 15 is half a WALLABY beam.
    DefaultRadSep="15.0"
    #   Set the default number of radial bins -- this is large to deal with some of the larger detections.
    DefaultNRad="50"
    #   Set the automatic fitting method.  For flat disks, a polynmial of degree 0 is preferable, but this is less important when the 2-stage switch is false.
    FitPolynDegree="0"
    #   Set the parameters that should be fit from the initial estimate by Barolo
    ParamsToFit="VROT\t INC\t PA\t  XPOS\t YPOS \t VSYS"
    #   Place all fitting options into a dictionary
    FittingOptions={'POLYN':FitPolynDegree,'FLAGERRORS':Errors,'TWOSTAGE':TwoStage,'PLOTS':PlotSwitch,'NORM':Normalization,'NRADII':DefaultNRad,'RADSEP':DefaultRadSep,'FREE':ParamsToFit,'VDISP':DefaultDispersion}
    return FittingOptions
    
def SetBaroloMaskingOptions(ObjDict,SoFiASwitch):
    """
        This function sets the masking options that 3D-Barolo will use
    """
    #   If the SoFiASwitch is false, use the built-in mask maker.
    if SoFiASwitch == False:
        MaskStr="\nMASK\t SEARCH\n"
        MaskStr+="SNRCUT\t 3.0\n"
        MaskStr+="GROWTHCUT\t 2.5\n"
    else:
    #   Otherwise use the SoFiA mask file -- note this will only work when using the full resolution WALLABY cubes
        MaskFile=ObjDict['ObjFileBaseName']+"_mask.fits.gz"
        MaskStr="\nMASK\t file("+MaskFile+")\n"

    return MaskStr
