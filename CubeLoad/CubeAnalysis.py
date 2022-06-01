import numpy as np

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

"""
    This module contains routines for loading and doing some basic cube analysis.  The routines contained in this module are:
    CubeAnalysisFuncDict --> This places particularly useful functions into a dictionary for access in other modules without explicit imports.
    SoFiA_CubeAnalysis --> This function loads and analyzes a specific galaxy's SoFiA cubelet detection.
    BasicCubeLoad --> This function loads in a cubelet.
    BasicCubeAnalysis --> This function loads in a cubelet and does a few basic analysis steps
    GetStartPixelLoc --> This function specifies the pixel at the start of a specific cubelet for later MCG analysis.
    TotalCubeHI -->  This function calculates the total HI in a cubelet
    MassCalcFromCube --> This function calculates the mass in a cubelet based on the total HI
    GetVelsFromHeader --> This function makes an array of velocities corresponding to each channel.
    ConstructPVDiagram --> This function calculates a basic PV diagram.
    ConstructModelBasedPVDiagram --> This function calculates a PV diagram based on a model cubelet.
    MakeQuickCubeMeasures --> This function calculates a number of useful measurements from a cube that are used in the diagnostic plotting functions
"""

def CubeAnalysisFuncDict():
    """
        This function stores other functions in this module into a dictionary.
    """
    CubeFuncs={'SoFiA_CubeAnalysis':SoFiA_CubeAnalysis,'TotalCubeHI':TotalCubeHI,'ConstructPVDiagram':ConstructPVDiagram,'MassCalcFromCube':MassCalcFromCube,'GetVelsFromHeader':GetVelsFromHeader,'GetStartPixelLoc':GetStartPixelLoc,'BasicCubeLoad':BasicCubeLoad,'BasicCubeAnalysis':BasicCubeAnalysis,'ConstructModelBasedPVDiagram':ConstructModelBasedPVDiagram,'MakeQuickCubeMeasures':MakeQuickCubeMeasures}
    
    return CubeFuncs


def SoFiA_CubeAnalysis(ObjDict):
    """
        This function loads in a SoFiA cube, mask, gets the mass, total HI, and constructs a major axis PV diagram.  All these quantities are saved to a dictionary for later tests and comparisons with models.
    """
    #   Set the names of the velocity cube and the mask cube
    VelCubeName=ObjDict['VelocityCubeFileName']
    MaskCubeName=ObjDict['MaskFileName']
    #   Load in the velocity cube
    CubeHeader,CubeData,CubeWCS=BasicCubeLoad(VelCubeName)
    #   Load in the Masked cube
    MaskCubeHeader,MaskData,MaskCubeWCS=BasicCubeLoad(MaskCubeName)
    #   Make the masked data cube
    MaskedCubeData=CubeData*MaskData
    #   Get the start location of the pixels from the header -- this is mainly for MCG comparisons
    CubeHeader=GetStartPixelLoc(CubeWCS,CubeHeader)
    
    #   Set the beam size
    BeamSize=CubeHeader['BMAJ']
    #   Set the pixel size
    PixSize=CubeHeader['CDELT2']
    #   Set the channel size in km/s
    ChannelSize=CubeHeader['CDELT3']
    ChannelSize=ChannelSize/1000.
    #   Get the beam size in units of pixels
    BeamSize_Pix=BeamSize/PixSize
    #   Get the distance from the measurement catalogue
    Distance=ObjDict['MeasureCat']['DISTANCE']
    #   Get the integrated HI of the masked cube
    IntegratedHI=TotalCubeHI(MaskedCubeData,PixSize,BeamSize)
    #   Get the mass from the cube
    CubeMass=MassCalcFromCube(MaskedCubeData,Distance,ChannelSize,BeamSize,PixSize)
    #   Get the velocities from the cube header
    CubeVels=np.array(GetVelsFromHeader(CubeHeader))
    
    #   Construct the major axis PV diagram
    #       This first requires setting the pixel center to the SoFiA center
    CenterPos=[[ObjDict['CatEntry']['ra'],ObjDict['CatEntry']['dec'],0]]
    CentPix=CubeWCS.wcs_world2pix(CenterPos,0)
    CentPix=[CentPix[0,0],CentPix[0,1]]
    #   The SoFiA PA needs a 180 degree adjustment to be consistent with kinematic measures.
    PA=ObjDict['CatEntry']['kin_pa']+180.
    #       Now the major axis PV can be calculated
    MajorAxisPV=ConstructPVDiagram(MaskedCubeData,PA,CentPix,BeamSize_Pix,CubeVels)
    #   Save the RMS value for tests later on.
    RMS=ObjDict['CatEntry']['rms']
    #   Place all the info into a dictionary and return that
    CubeInfo={'Data':CubeData,'Mask':MaskData,'MaskedCube':MaskedCubeData\
        ,'CubeHeader':CubeHeader, 'Distance':Distance,'IntegratedHI':IntegratedHI\
        ,'Mass':CubeMass,'CubeWCS':CubeWCS,'Noise':RMS,'CubeVels':CubeVels,'PV':MajorAxisPV,'CentPix':CentPix
        }
    return CubeInfo
    
    
def BasicCubeLoad(CubeName):
    """
        This function uses the astropy fits library to load in a cube and seperate it into a header, data, and WCS object.
    """
    #   Open up the cube
    CubeHDU=fits.open(CubeName)
    #   Load in the header
    CubeHeader=CubeHDU[0].header
    #   Load in the data
    CubeData=CubeHDU[0].data
    #   Get the wcs values
    CubeWCS = wcs.WCS(CubeHeader)
    #   Close the cube
    CubeHDU.close()
    return CubeHeader,CubeData,CubeWCS
    
def BasicCubeAnalysis(CubeName):
    """
        This function loads in the cube information using BasicCubeLoad, calculates the velocity channel array and starting pixel locations.  All this is stored into a dictionary.  It can be used on any particular cube provided the file name is specified.  This is different than the SoFiA analysis, which requires a full galaxy object dictionary to run.
    """
    #   Load in the cube header
    CubeHeader,CubeData,CubeWCS=BasicCubeLoad(CubeName)
    #   Get the cube velocities
    CubeVels=np.array(GetVelsFromHeader(CubeHeader))
    #   Get the start location of the pixels from the header
    CubeHeader=GetStartPixelLoc(CubeWCS,CubeHeader)
    #   Store all the info into a dictionary and return in
    CubeInfo={'Data':CubeData,'CubeHeader':CubeHeader,'CubeWCS':CubeWCS,'CubeVels':CubeVels}
    return CubeInfo
    
def GetStartPixelLoc(CubeWCS,CubeHeader):
    """
        This function gets the RA and DEC coordinates for the first pixel in a cube.  It used for some MCG tests but is not really required for the full WKAPP analysis.
    """
    #   Set the pixel coordinate as 0,0,0
    pixcrd=[[0,0,0]]
    #   Get the Equatorial coordinates of the pixel
    RealCoord=CubeWCS.wcs_pix2world(pixcrd,0)
    #   Store the spatial coordinates in the header.
    CubeHeader['CRLOC1']=RealCoord[0,0]
    CubeHeader['CRLOC2']=RealCoord[0,1]
    return CubeHeader

def TotalCubeHI(cube,PixSize,BeamSize):
    """
        This function calculates to total HI in a cube  in terms of Jy.
    """
    cubeTot=np.nansum(cube)         #Get the total signal in cube units (Jy/beam)
    BeamArea=np.pi*BeamSize**2./(4.*np.log(2.)) #Get the beam area in degrees
    pixperbeam=BeamArea/PixSize**2. #Find the number of pixels per beam
    HISignal=cubeTot/pixperbeam #Get the integrated HI signal in Jy
    return HISignal


def MassCalcFromCube(cube,Distance,ChannelSize,BeamSize,PixSize):
    """
        This function calculates the HI mass in the cube from the flux based on ()
    """
    #   Get the total signal in Jy/beam
    Signal_JyBeam=np.nansum(cube)
    #       Get the beam size in degrees
    beamarea=(np.pi*BeamSize**2.)/(4.*np.log(2.))
    #   Ge the number of pixels/beam
    pixperbeam=beamarea/(abs(PixSize)*abs(PixSize))
    #   Get the signal in Jy
    Signal_Jy=Signal_JyBeam/pixperbeam
    #   Get the mass
    Mass=0.236*(Distance*1000.)**2.*Signal_Jy*abs(ChannelSize)  #Original version
    return Mass


def GetVelsFromHeader(CubeHeader):
    """
        This function gets the velocity for each channel in a cube.  This is used for calculations of moment maps and PV diagrams.
    """
    #   Make an empty list of velocities.
    Vels=[]
    #   Loop through all channels
    for i in range(CubeHeader['NAXIS3']):
        #   Get the velocity at a specific channel
        VTemp=float(i-CubeHeader['CRPIX3']+1)*CubeHeader['CDELT3']+CubeHeader['CRVAL3']
        #   Append the channel velocity to the velocity list.
        Vels.append(VTemp)
    return Vels


def ConstructPVDiagram(CubeData,Angle,PixelCenter,BeamSize,Vels):
    """
        This function calculates a PV diagram for a cube along a specified axis angle.  It does a reasonable quick calculation, but there are better versions that are symmetric about the pixel center, which this version is not.
    """
    #   First adjust the angle from the normal PA reference frame (from the +Y axis) to be from the X-axis.  Also put the angle into radians
    #   Adjust the angle to be measured from the +X axis instead of the +Y axis.
    AngUse=Angle+90.
    #   Make sure the angle is between 0 and 2 pi
    if AngUse>360.:
        AngUse=AngUse-360.
    if AngUse<0.:
        AngUse=AngUse+360
    AngUse=AngUse*np.pi/180.
    #   Set up an empty array based on the cube.
    PV=np.zeros([np.shape(CubeData)[2],np.shape(CubeData)[0]])
    
    #   Loop through all spatial pixels in the cube.
    for i in range(np.shape(CubeData)[1]):
        for j in range(np.shape(CubeData)[2]):
            #   Turn the pixel coordinates into relative coordinates from the pixel center.
            X=(j-PixelCenter[0])
            Y=(i-PixelCenter[1])
            #   Rotate the pixels to be along the positive x-axis.
            XP=X*np.cos(-AngUse)-Y*np.sin(-AngUse)
            YP=X*np.sin(-AngUse)+Y*np.cos(-AngUse)
            #   Get the indices of the rotated pixels
            k=int(round(XP+PixelCenter[0]))
            l=int(round(YP+PixelCenter[1]))
            #   Check that the pixels along the specified angle are inside the array.
            if k >= 0 and k < np.shape(CubeData)[2]:
                #   Only include pixels within half a beam of the specified axis (note that the beam size must be in pixels).
                if abs(YP) <= BeamSize/2.:
                    #   Loop through the velocity channels
                    for m in range(np.shape(CubeData)[0]):
                        #   Add in the flux from cells to the PV diagram 'pixels'
                        if np.isnan(CubeData[m,i,j]):
                            PV[k,m]+=0.
                        else:
                            PV[k,m]+=CubeData[m,i,j]
    #   Return the PV cube array.
    return PV
    
def ConstructModelBasedPVDiagram(CubeData,Angle,AvgModel,BeamSize,OriDataCube):
    """
        This function calculates the PV diagram for a specific model
    """
    #   copy in the velocity array from the original cube.
    Vels=OriDataCube['CubeVels']
    #   Adjust the angle to be measured from the +X axis instead of the +Y axis.
    AngUse=Angle+90.
    #   Make sure the angle is between 0 and 2 pi
    if AngUse>360.:
        AngUse=AngUse-360.
    if AngUse<0.:
        AngUse=AngUse+360
    AngUse=AngUse*np.pi/180.

    #   Get the size of each bin
    dR=AvgModel['R_SD'][1]-AvgModel['R_SD'][0]
    #   Get the largest extent possible from the model
    RTest=AvgModel['R_SD'][-1]+2*dR
    #   Get the number of radial pixels that covers the extent of the model
    nRPix=int(RTest/np.abs(OriDataCube['CubeHeader']['CDELT1']*3600.))
    #   Set the number of spatial pixels required to cover the full size of the model
    nSpatialPix=2*nRPix+1
    #   Allocate the PV diagram using the number of spatial pixels required to place the center in the center of the array and the size of the velocity dimension
    PV=np.zeros([nSpatialPix,np.shape(CubeData)[0]])
    #   Set up a pair of 'center' pixel coordinates.
    PixelCenter=np.zeros(2)
    PixelCenter[0]=AvgModel['XCENTER'][0]
    PixelCenter[1]=AvgModel['YCENTER'][0]
    
    #   Loop through all spatial pixels in the cube.
    for i in range(np.shape(CubeData)[1]):
        for j in range(np.shape(CubeData)[2]):
            #   Turn the pixel coordinates into relative coordinates from the pixel center.
            X=(j-PixelCenter[0])
            Y=(i-PixelCenter[1])
            #   Rotate the pixels to be along the positive x-axis.
            XP=X*np.cos(-AngUse)-Y*np.sin(-AngUse)
            YP=X*np.sin(-AngUse)+Y*np.cos(-AngUse)
            #   Get the indices of the rotated pixels
            k=int(round(XP+nSpatialPix/2))
            l=int(round(YP+nSpatialPix/2))
            #   Check that the pixels along the specified angle are inside the array.
            if k >= 0 and k < nSpatialPix:
                #   Only include pixels within half a beam of the specified axis (note that the beam size must be in pixels).
                if abs(YP) <= BeamSize/2.:
                    #   Loop through the velocity channels
                    for m in range(np.shape(CubeData)[0]):
                        #   Add in the flux from cells to the PV diagram 'pixels'
                        if np.isnan(CubeData[m,i,j]):
                            PV[k,m]+=0.
                        else:
                            PV[k,m]+=CubeData[m,i,j]
    #   Return the PV array
    return PV



def MakeQuickCubeMeasures(ObjDict,CubeInfo,AnalysisFncs):
    """
        This function makes a set of quick measurements that are added to the plot in the various labels.
    """
    #   Start by making an empty dictionary
    CubeMeasureDict={}
    #   Get the diameter of the object in beams from the measurement catalogue (see Wallaby_BasicAnalysis for all the parts of this dictionary.
    CubeMeasureDict['Diameter']=ObjDict['MeasureCat']['SIZE'].values[0]
    #   Set the radius in and diameter in arcseconds rather than beams.
    CubeMeasureDict['Radius']=CubeMeasureDict['Diameter']/2.*30.
    CubeMeasureDict['DiameterR']=CubeMeasureDict['Diameter']*30.
    #   Use the read in mass from the cube (see CubeLoad/CubeAnalysis.py)
    CubeMeasureDict['logM']=np.log10(CubeInfo['Mass'].values[0])
    CubeMeasureDict['Mass']=CubeInfo['Mass'].values[0]
    #   Set the distance from the measurement catalogue.
    CubeMeasureDict['Distance']=ObjDict['MeasureCat']['DISTANCE'].values[0]
    #   Calculate the predicted value for R_HI from the HI mass-size relation
    RHITest=AnalysisFncs['AstroFncs']['CalcRHI'](CubeMeasureDict['logM'])
    #   Convert this size to arcseconds for use in the diagnostic plots.
    CubeMeasureDict['RHI']=AnalysisFncs['AstroFncs']['ProjectedSize_ArcSecond'](RHITest,CubeMeasureDict['Distance'])
    #   Store the central pixel from the cube
    CubeMeasureDict['CentPix']=CubeInfo['CentPix']

    return CubeMeasureDict
