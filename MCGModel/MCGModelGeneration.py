import numpy as np
import astropy
import copy as copy
import os

"""
    This module contains routines needed to convert a tilted ring dictionary into a model cube using the MockCubeGenerator (MCG) program from MCGSuite.  It contains the routines:
    MCGModelFnc --> This function stores some of the routines into a dictionary for easy passing to different routines.
    MCG_Defs --> This function contains some definitions for running MCG
    MakeMCGModel --> This function runs MCG to make a particular tilted ring model
    MakeMCG_MainFile --> This function writes the main input file for MCG
    MakeMCGDataCubeHeader --> This function writes the file that MCG uses to determine the target datacube parameters
    MakeMCG_TiltedRingInput --> This function writes the tilted ring model file for MCG
"""

def MCGModelFnc():
    """
        This function stores some of the routines into a dictionary for easy passing to different routines.
    """
    MCGFncDict={'MakeMCGModel':MakeMCGModel,'MCG_Defs':MCG_Defs}
    return MCGFncDict


def MCG_Defs(ObjName):
    """
        This function contains some definitions for running MCG
    """
    #   First hard code in the location of the MockCubeGenerator program (MCG)
    MCGPath="/Users/nate/Dropbox/GitProjects/21-cm-Emission/PSOFT14/MCGSuite/Programs/MockCubeGenerator"
    #   Set the names of the three MCG input files
    MCG_DatacubeHeaderName="TempMCG_DCHeader.in"
    MCG_TiltedRingFileName="TempMCG_TR.in"
    MCG_MainInFile="TempMCG.in"
    #   Name the model after the galaxy
    MCG_ModelName=ObjName
    #   Select the noiseless cube name to be used for later comparisons
    FinalModelCube=MCG_ModelName+"/"+MCG_ModelName+"_ConvolvedSourceCube.fits"
    #   Store all the information into a dictionary for ease of access
    MCGNames={'MCGPath':MCGPath ,'OutputName':MCG_ModelName,'MainFile':MCG_MainInFile\
        ,'DatacubeFile':MCG_DatacubeHeaderName,'TiltedRingFile':MCG_TiltedRingFileName,'CubeFile':FinalModelCube}
    return MCGNames


def MakeMCGModel(MCGDict,ModelParams,CubeHeader):
    """
        This function makes a mock cube for a given tilted ring model (ModelParams) using MCG.
    """
    #       This function requires a set of Tilted ring model parameters and a fits cube header
    #   Make the MCG Main file
    MakeMCG_MainFile(MCGDict)
    #   Make the MCG Data cube header file
    MakeMCGDataCubeHeader(MCGDict['DatacubeFile'],CubeHeader)
    #   Now make the tilted ring input file
    MakeMCG_TiltedRingInput(MCGDict['TiltedRingFile'],ModelParams)
    #   Now generate the MCG Model
    MakeCubeCmd=MCGDict['MCGPath']+" " + MCGDict['MainFile']
    os.system(MakeCubeCmd)
    #   Once done, remove the input files
    CleanCmd="rm " + MCGDict['MainFile'] +" " + MCGDict['DatacubeFile'] + " " + MCGDict['TiltedRingFile']
    os.system(CleanCmd)

def MakeMCG_MainFile(MCGNames):
    """
        This function writes the main input file needed for MCG.  Since MCG is written in Fortran, the input file format is fixed and not particularly flexible.
    """
    #   Open up the input file for writting
    file=open(MCGNames['MainFile'],"w")
    #   Write out the name of the output folder (named after the galaxy)
    file.write("#    Base name for the output folder and all file names \n")
    file.write(MCGNames['OutputName']+"\n")
    #   There are different output options for MCG...2 gives the most outputs
    file.write("#    Minimal, Moderate, of Full Outputs (0,1, and 2 respectively)\n")
    file.write("2 \n")
    #   Give the name of the MCG datacube definitions input file
    file.write("#    Name of the file containing the cube and beam definitions\n")
    file.write(MCGNames['DatacubeFile']+"\n")
    #   Give the name of the MCG Tilted Ring model input file
    file.write("#    Name of the file containing the underlying ring model\n")
    file.write(MCGNames['TiltedRingFile']+"\n")
    #   Set the units for the noise
    file.write("#    Noise Unit Switch (0=mJy/beam)\n")
    file.write("0 \n")
    #   Set the noise value
    file.write("#    Noise Value\n")
    file.write("1.6 \n")
    #   Set the random seed -- any positive number uses the system time
    file.write("#    Random Seed\n")
    file.write("1 \n")
    #   Close the main input file.
    file.close()

def MakeMCGDataCubeHeader(MCGFileName,CubeHeader):
    """
        This function writes the MCG file that determines the cube parameters.  It uses the existing cube header so the mock data cube has the same properties/resolution as the observed cube.
    """
    #   First open the file
    file=open(MCGFileName,"w")
    #   Next write out the cube dimensions
    file.write("#    number of pixels and channel\n")
    file.write(str(CubeHeader['NAXIS1'])+"\t" +str(CubeHeader['NAXIS2']) \
               +"\t" +str(CubeHeader['NAXIS3']) +"\n")
    #   Then add in the units that MCG will understand
    file.write("#    Pixel Size Units (0=degrees, 1=arcsec)  - Channel Size Units (0=m/s, 1=km/s) - Reference Pixel Units (0=degrees,1=arcsec) - Reference Channel Units (0=m/s, 1=km/s)  — Beam axis Units (0=arcsec) — Beam Rotation Angle Units (0=degrees)\n")
    file.write("0\t 0 \t 0 \t 0 \t 1 \t 0\n")
    #   Now write the size of the pixels/channels
    file.write("#    Pixel dimensions and channel size\n")
    file.write(str(CubeHeader['CDELT1']) +"\t" +str(CubeHeader['CDELT2'])\
               +"\t"+str(CubeHeader['CDELT3']) +"\t"+" \n")
    #   And give the reference location for the pixels
    file.write("#    Reference Pixel(channel) in each dimension (CRPIX values)\n")
    file.write("1" +"\t" +"1"\
               +"\t"+str(CubeHeader['CRPIX3']) +"\t"+" \n")
    #   Set the value in RA, DEC, and m/s for each dimension
    file.write("#    Reference value in each dimension (CRVAL values)\n")
    file.write(str(CubeHeader['CRLOC1']) +"\t" +str(CubeHeader['CRLOC2'])\
               +"\t"+str(CubeHeader['CRVAL3']) +"\t"+" \n")
    #   Write out the beam parameters for smoothing
    file.write("#    Beam dimensions (major, minor, position angle)\n")
    file.write(str(CubeHeader['BMAJ']*3600.) +"\t" +str(CubeHeader['BMIN']*3600.)
               +"\t"+str(CubeHeader['BPA']) +"\t"+" \n")
    #   The beam can be smoothed over some number of sigma before truncation...5 numerically includes all the flux.
    file.write("#    number of sigma lengths to reach with the beam\n")
    file.write("5. \n")
    #   It is possible to include instrumental velocity smoothing, but this isn't necessary for WALLABY fits.
    file.write("#    Type of velocity smoothing to use (0=none, 1=Gaussian)\n")
    file.write("0 \n")
    #   If we were smoothing, set the velocity sigma.  Due to the fixed nature of MCG, something needs to exist here.
    file.write("#    The velocity smoothing sigma if using Gaussian smoothing (km/s)\n")
    file.write("10. \n")
    #   Close the file.
    file.close()

def MakeMCG_TiltedRingInput(MCGFileName,ParamDict):
    """
        This function writes the tilted ring model input file that MCG uses to build the model cube.
    """
    #   Set the total number of rings in the model -- note that for MCG, the SD and VRot profiles must have the same extent.
    nRings=np.shape(ParamDict['R'])[0]
    #   Set the ring width -- if there's more than one ring, use the difference between the 2nd and 1st ring
    try:
        Rwidth=ParamDict['R'][1]-ParamDict['R'][0]
        #   Otherwise, just double the ring radius (as MCG uses it as a midpoint)
    except:
        Rwidth=2.*ParamDict['R'][0]
    #   Open up the input file
    file=open(MCGFileName,"w")
    #   Write out the number of rings in the model
    file.write("#   Number of Rings in the model\n")
    file.write(str(nRings)+"\n")
    #   Write out the way the particle density and flux is set -- 0 is the default
    file.write("#    The cloud mode you are using\n")
    file.write("0 \n")
    #   Get the number of particles per area
    file.write("#    The base cloud surface density\n")
    file.write("100. \n")
    #   Set the units for different parameters
    file.write("#    Central Position Switch (0=degrees, 1=arcsec,2=pixels), Inc/PA Unit switch (0=degrees, 1=arcsec), Velocity Units (0=m/s, 1=km/s), Brightness Units (0=Jy km/s arcsec^-2)\n")
    file.write("2 \t 0 \t 1 \t  0 \n")
    #   Set the model header
    file.write("#    The parameters in each radial bin\n")
    file.write("#    Rmid    Rwidth    Xcent    Ycent    Inc    PA    VSys    VRot    VRad    Vvert    VDisp    dvdz    Sigma        z0    zGradStart\n")
    #    Loop through all rings
    for i in range(nRings):
        #   Write the ring parameters needed for MCG.  For certain parameters that aren't modelled like radial velocities, set them to zero
        file.write(str(ParamDict['R'][i])+"\t"+str(Rwidth) +"\t" + str(ParamDict['XCENTER'][0]) \
                   + "\t"+ str(ParamDict['YCENTER'][0]) +"\t"\
                   + str(ParamDict['INCLINATION'][i])\
                   + "\t" +str(ParamDict['POSITIONANGLE'][i]) +"\t"\
                   + str(ParamDict['VSYS'][i])+ "\t" + str(ParamDict['VROT'][i]) \
                   +"\t" + str(ParamDict['VRAD'][i]) + "\t" + "0." +"\t" \
                   + str(ParamDict['VDISPERSION'][i]) + "\t"+ "0." + "\t" \
                   + str(ParamDict['SURFDENS_FACEON'][i])+"\t" + str(ParamDict['Z0'][i])\
                   + "\t" + str(5.*ParamDict['Z0'][i]) +   "\n")
    #   Close the file
    file.close()
