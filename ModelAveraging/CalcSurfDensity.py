import numpy as np
import astropy
import copy as copy
import os
from scipy import interpolate
"""
    This module contains the routines that generate the average surface density from a moment map average model geometric parameters.  It contains the routines:

    GetNewSurfDens --> This function constructs the surface brightness for the new average model using the average geometry
    CalculateEllipticalRadius --> This function calculates the elliptical radius of a point
    DetermineREllipIndx --> This function gets the index for an elliptical radius for a point
    DetermineRMax --> This function determines the largest radius the surface density may reach
    GetPixelPoint --> This function converts a pixel into a point
    ConstructMom0 --> This function constructs a moment 0 map for the masked cube to be used for the surface density
"""

def GetNewSurfDens(AvgModel,CubeInfo,AnalysisFncs):
    """
        This function calculates the surface density for an tilted ring model from the moment 0 map.
    """
    #   First get the moment 0 map from the data using the masked cube
    #       Construct a moment 0 map
    Mom0=ConstructMom0(CubeInfo,AnalysisFncs)

    #   Get the ellipticity from the inclination
        #   First get the inclination and position angle in radians
    IncRad=AvgModel['INCLINATION'][0]*np.pi/180.
    PARad=(AvgModel['POSITIONANGLE'][0]+90.)*np.pi/180.
        #   Now get the ellipticity
    Ellipticity=np.cos(IncRad)
    #   Next we need to get the size of the surface density radial grid using the DetermineRMax function
    maxR,dR=DetermineRMax(AvgModel,CubeInfo,PARad,Mom0)
    #   It's necessary to have the radial grid size in pixels rather than arcseconds for some calculations.
    dRPix=dR/(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)
    #   Initialize the summation arrays needed for the SD calculation
    FluxTot=np.zeros(maxR)
    nPixTot=np.zeros(maxR)
    FluxDiff=np.zeros(maxR)
    AvgModel['R_SD']=np.zeros(maxR)
    #   Set up the radial grid in terms of arcseconds
    for i in range(maxR):
        AvgModel['R_SD'][i]=AvgModel['R'][0]+i*dR
    #   Set the first midpoint between the inner circle and the first ring in pixels as this is necessary for indexing
    RMid1=(AvgModel['R'][0]+AvgModel['R'][1])/2.
    #   And, as usual, this requires a value in pixels
    RMid1Pix=RMid1/(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)
    #       Loop through the full moment0 map
    for i in range(np.shape(Mom0)[0]):
        for j in range(np.shape(Mom0)[1]):
            #       Get the elliptical radial index for the current pixel
            RadIndx=DetermineREllipIndx(i,j,AvgModel,PARad,Ellipticity,RMid1Pix,dRPix)
            #   Check if the radius is inside the maximum size
            if RadIndx < maxR:
                #   Check if the Mom0 is a NaN
                if np.isnan(Mom0[i,j]):
                    FluxTot[RadIndx]+=0.
                    #   If it's not, add the moment 0 flux to the calculation
                else:
                    FluxTot[RadIndx]+=Mom0[i,j]
                nPixTot[RadIndx]+=1
    #   The pixel areay is needed for the SD calculation
    PixArea=(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)**2.
    #   Divide the flux by the area to get the SD
    AvgModel['SURFDENS']=copy.copy(FluxTot/(nPixTot*PixArea))#*np.cos(IncRad)
    #   Now get the uncertainty
    #       Again, loop through all pixels
    FluxAvg=copy.copy(FluxTot/(nPixTot))
    for i in range(np.shape(Mom0)[0]):
        for j in range(np.shape(Mom0)[1]):
            #       Get the elliptical radial index for the current pixel
            RadIndx=DetermineREllipIndx(i,j,AvgModel,PARad,Ellipticity,RMid1Pix,dRPix)
            #   Check if the radius is inside the maximum size
            if RadIndx < maxR:
                #   Check if the Mom0 is a NaN
                if np.isnan(Mom0[i,j]):
                    FluxDiff[RadIndx]+=(0.-FluxAvg[RadIndx])**2.
                    #   If it's not, add the square of the difference between the pixel flux and average radial flux
                else:
                    FluxDiff[RadIndx]+=(Mom0[i,j]-FluxAvg[RadIndx])**2.
    #   Calculate the surface density uncertainty using the standard deviation
    #       First get the number of beams in the ring -- since dRPix is half a beam, get the beam area via 2*pi*(FWHM/2.355)**2.
    BeamArea_Pix=(2.*np.pi*(dRPix*2/2.355)**2.)
    #       Now get the number of beams in the ring
    nBeams=nPixTot/BeamArea_Pix
    #       Make sure that the number of beams per ring is at least 1
    for i in  range(len(nBeams)):
        if nBeams[i] <=1:
            nBeams[i]=1
    #       Get the st.dev. as the variance/sqrt(number of beams)
    AvgModel['SURFDENS_ERR']=copy.copy(np.sqrt(FluxDiff/(nPixTot))/PixArea)/np.sqrt(nBeams)
    #   Convert the SD and uncertainty into Msol/pc^2 units
    SDConv=1.24756e+20/(6.0574E5*1.823E18*(2.*np.pi/np.log(256.)))
    AvgModel['SDCONV']=SDConv
    AvgModel['SURFDENS']/=SDConv
    AvgModel['SURFDENS_ERR']/=SDConv
    #   Do a face on correction if useful
    AvgModel['SURFDENS_FACEON']=AvgModel['SURFDENS']*np.cos(IncRad)
    AvgModel['SURFDENS_FACEON_ERR']=AvgModel['SURFDENS_ERR']*np.cos(IncRad)
    return AvgModel


def DetermineREllipIndx(i,j,AvgModel,PARad,Ellipticity,RMid1Pix,dRPix):
    """
        This function determines the index of a point in elliptical coordinates on the surface density radial grid.
    """
    #   Turn the pixel coordinates into X and Y
    Y=copy.copy(i)
    X=copy.copy(j)
    #   Get the elliptical radius for this point
    REllip=CalculateEllipticalRadius(X,Y,AvgModel,PARad,Ellipticity)
    #   Get the index for this bin using the size of the bins in pixels
    RealRadBin=(REllip-RMid1Pix)/dRPix
    #   Check whether we managed to get an index less than zero
    if RealRadBin <= 0.:
        RadIndx=0
    else:
        #   Convert the real into an integer
        RadIndx=int((REllip-RMid1Pix)/dRPix)+1
    #   Return the radial index
    return RadIndx
    
def CalculateEllipticalRadius(X,Y,AvgModel,PARad,Ellipticity):
    """
        This function gets the elliptical radius of a point given a specific position angle and ellipticity
    """
    #   First get the pixel position relative to the model center
    XCent=X-AvgModel['XCENTER'][0]
    YCent=Y-AvgModel['YCENTER'][0]
    #   Next rotate the points using the PA in radians
    XRot=XCent*np.cos(-PARad)-YCent*np.sin(-PARad)
    YRot=XCent*np.sin(-PARad)+YCent*np.cos(-PARad)
    #   Get the Yprime coorddinate by divinding by the ellipticity
    YEllip=YRot/Ellipticity
    #   Calculate the elliptical radius
    REllip=np.sqrt(XRot**2.+YEllip**2.)
    return REllip

def DetermineRMax(AvgModel,CubeInfo,PARad,Mom0):
    """
        This function gets the maximum radius for the surface density profile
    """
        #   Get the radial size of the bins in pixels
    dR=AvgModel['R'][1]-AvgModel['R'][0]
    #       Set the maximum radius by setting it to the length of the diaganol across the image
    maxX=(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)*CubeInfo['CubeHeader']['NAXIS1']
    maxY=(abs(CubeInfo['CubeHeader']['CDELT2'])*3600.)*CubeInfo['CubeHeader']['NAXIS2']
    maxR=int(np.sqrt(maxX**2.+maxY**2)/dR)
    
    #   Loop through the maximum possible length of the grid
    for i in range(maxR):
        #   Get the radius of the current point
        TestR=(AvgModel['R'][0]+(i+0.5)*dR)
        #   Convert it to a point in pixels
        TestR=TestR/(abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)
        #   Turn it to a point in X and Y for the pixels
        X,Y=GetPixelPoint(TestR,AvgModel,PARad)
        #   Only start checking the limits once we reach the end of the existing rotation curve grid
        if i >= len(AvgModel['R'])-1:
            #   Check whether X and Y are inside the image
            if X <0 or X >CubeInfo['CubeHeader']['NAXIS1']-1:
                print("Beyond X Lims SD break")
                break
            if Y <0 or Y >CubeInfo['CubeHeader']['NAXIS2']-1:
                print("Beyond Y Lims SD break")
                break
            #   If it's inside, get the indices of the X and Y values
            iTest=int(Y)
            jTest=int(X)
            #   Check whether the map has been masked out
            if Mom0[iTest,jTest] <=0.:
                print("Flux Lims SD break", X,Y,Mom0[iTest,jTest])
                break
    #   Once the loop ends (either from hitting the masked out points or the edge of the profile, i+1 will be the maximum number of radial points to use in the SD profile calculation.
    maxR=i+1
    #   Return both the maximum number of grid points and the size of the grid.
    return maxR,dR
    
    
def GetPixelPoint(R,AvgModel,PARad):
    """
        This function turn a radius + position angle into a pixel location
    """
    #   In elliptical coordinates, the point is at (R,0)
    XEllip=R
    YEllip=0.
    #   Rotate this to physical coordinates
    XRot=XEllip*np.cos(PARad)-YEllip*np.sin(PARad)
    YRot=XEllip*np.sin(PARad)+YEllip*np.cos(PARad)
    #   Go from coordinates relative to the center to absolut pixel coordinates.
    XPos=XRot+AvgModel['XCENTER'][0]
    YPos=YRot+AvgModel['YCENTER'][0]
    #   Return the position.
    return XPos,YPos

def ConstructMom0(Cube,AnalysisFncs):
    """
        This function makes the moment 0 map for the cube.  It uses the moment map function from DiagnosticPlots/MomentMapFncs.py
    """
    #   Set the 'MaskedData' keyword to the 'MaskedCube' keyword -- the moment map calculation using MaskedData instead of MaskedCube, but the two are the same.
    Cube['MaskedData']=Cube['MaskedCube']
    #   Make the moment 0 map
    Mom0=AnalysisFncs['PlotFncs']['MakeMomMap'](Cube,0,1)
    #   The moment 0 map needs to be converted to the correct units using the cube header
    CubeHeader=Cube['CubeHeader']
    #   Convert to Jy/beam km/s
    Mom0=(Mom0*abs(CubeHeader['CDELT3']/1000.))
    #   Convert to Jy/pixel km/s
    BeamSize=CubeHeader['BMAJ']/abs(CubeHeader['CDELT1'])
    BeamArea3=2.*np.pi*(BeamSize/2.355)**2.
    Mom0=Mom0/BeamArea3
    #   Return the moment 0 map
    return Mom0
