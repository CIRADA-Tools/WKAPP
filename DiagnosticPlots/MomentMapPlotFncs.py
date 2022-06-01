import numpy as np
import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs
import copy as copy

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator

from matplotlib.patches import Ellipse

"""
    This module contains a number of plotting functions that are useful for both diagnostic and average model plots.  This module contains routines focused on making moment maps.  It contains the routines:

    MomentPlot --> This function makes a basic moment map panel for a plot
    AddCenterToMomMap --> This function adds a point to a moment map panel indicating the center
    AddArrowToMomMap --> This function adds an arrow to a moment map panel indicating the position angle
    AddVelContoursToMomentPlot --> This function adds a set of velocity contours to a moment map based on some other moment map.
    MakeMomMap --> This function makes a set of moment maps for some cube.
    QuickMaskCube --> This does a quick masking of a data cube by a galaxy's SoFiA mask cube.
    
    DetermineExtent --> This function determines the extent for a moment map.
    FormatMomMap --> This function formats and cleans up the moment map panels.
"""

def QuickMaskCube(ObsCube,ObjDict):
    """
        This function masks the observed data cube by the galaxy's SoFiA mask file and stores it in the ObsCube dictionary.
    """
    #   First get rid of any nan's
    ObsCube['Data']=np.nan_to_num(ObsCube['Data'])
    #   Now open up the SoFiA masked cube --> note that this means the 'observed cube' must be at the native 4 km/s resolution.
    MaskCube=fits.open(ObjDict['MaskFileName'])
    #   Create a 'MaskedData' array for use in making moment maps
    ObsCube['MaskedData']=ObsCube['Data']*MaskCube[0].data
    #   Close the mask file
    MaskCube.close()
    #   Return the observed cube dictionary.
    return ObsCube

def MomentPlot(fig,placement,Moment,ObjDict,ObsCube,XC,YC,AvgModel=False):
    """
        This function makes a moment map panel on a larger canvas.
    """
    #   First set up the moment map panel on the canvas using the placement rectangle.
    ax=fig.add_axes(placement)
    #   Next mask the 'observed cube' by the SoFiA mask -- the observed cube must be at 4 km/s resolution.
    ObsCube=QuickMaskCube(ObsCube,ObjDict)
    #   Make the appropriate moment map (indicated by the Moment integer).
    MomMapData=MakeMomMap(ObsCube,Moment,1)
    
    #   Set the base minimum and maximum values colormap on the map
    minval=np.nanmin(MomMapData)
    maxval=1.05*np.nanmax(MomMapData)
    #   Adjust the moment 1 colormap limits
    if Moment == 1:
        #   If we don't have a TR model to use when making the map, use the cube size
        minval=ObsCube['CubeVels'][-1]/1000.
        maxval=ObsCube['CubeVels'][0]/1000.
        #   If we do have a model, use the inclination corrected velocity to get the limits on the colormap
        if AvgModel != False:
            VS=AvgModel['VSYS'][0]
            dV=np.nanmax(AvgModel['VROT'])*np.sin(AvgModel['INCLINATION'][0]*np.pi/180.)
            minval=VS-1.1*dV
            maxval=VS+1.1*dV
    #   Set the colored maps
    if Moment == 0:
        Cmap='viridis'
    elif Moment == 1:
        Cmap='jet'
    elif Moment ==2 :
        Cmap='plasma'
    #   Format the moment map panel
    ExtentDict=FormatMomPlot(ax,Moment,XC,YC,ObsCube)
    #   Make the basic moment map using imshow to keep the proper aspect ratio
    ax.imshow(MomMapData,cmap=Cmap,vmin=minval,vmax=maxval,origin='lower',extent=ExtentDict['Extent'])
    #   For the moment 0 map, add in an ellipse to represent 1 beam
    if Moment ==0:
        #   First set the size of the beam from pixels to arcseconds
        ESize=5.*abs(ExtentDict['dX'])
        #   Set the center of the ellipse to pixel (7,7)
        Loc=7.
        #   Adjust the location to be coordinates in arcseconds using the 'ExtentDict' generated in DetermineExtent via FormatMomPlot
        LocX=(Loc-ExtentDict['RefPixX'])*ExtentDict['dX']+ExtentDict['RefValX']
        LocY=(Loc-ExtentDict['RefPixY'])*ExtentDict['dY']+ExtentDict['RefValY']
        #   Get the ellipse
        Ell=Ellipse([LocX,LocY], ESize, ESize, 0.,edgecolor='cyan',facecolor='none',lw=2)
        #   Add the ellipse to the map.
        ax.add_patch(Ell)

    return ax,MomMapData
    
def DetermineExtent(XC,YC,ObsCube):
    """
        This function sets the extent for moment maps in arcseconds
    """
    #   Set up the X and Y definitions needed
    #       First get the reference pixel
    RefPixX=XC
    #       Set the center to 0
    RefValX=0.
    #       Get the size of each pixel in arcseconds -- this assumes CDELT is in arcseconds
    dX=ObsCube['CubeHeader']['CDELT1']*3600.
    #       Get the number of pixels in the x direction
    nX=ObsCube['CubeHeader']['NAXIS1']
    #   Repeat the steps for the Y direction
    RefPixY=YC
    RefValY=0.
    dY=ObsCube['CubeHeader']['CDELT2']*3600.
    nY=ObsCube['CubeHeader']['NAXIS2']
    
    #  Get the lower and upper limits in arcseconds for both X and Y
    LowX=(-RefPixX)*dX+RefValX
    HighX=((nX-1)-RefPixX)*dX+RefValX
    LowY=(-RefPixY)*dY+RefValY
    HighY=((nY-1)-RefPixY)*dY+RefValY
    #   Set the extent array to be (X_low,X_high,Y_low,Y_high)
    Extent=[LowX,HighX,LowY,HighY]
    #   Store the extents and calculated values into a dictionary in case the extra info is needed elsewhere.
    ExtentDict={'Extent':Extent,'nX':nX,'RefPixX':RefPixX,'RefValX':RefValX,'dX':dX,'nY':nY,'RefPixY':RefPixY,'RefValY':RefValY,'dY':dY}
    #   Return the extent dictionary
    return ExtentDict
        
def FormatMomPlot(ax,Moment,XC,YC,ObsCube):
    """
        This function formats the moment map panel
    """
    #   Label the panel by the type of moment map it is showing
    PlotName="Mom"+str(Moment)
    #   Add the label to the canvas
    ax.text(0.5, 1.05,PlotName,fontsize=15,horizontalalignment='center',transform=ax.transAxes)
    #   Add minor axis ticks
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=4))
    #   Figure out how large the panel is in arcseconds from the center
    ExtentDict=DetermineExtent(XC,YC,ObsCube)
    #   Set the limits to the x and y coordinates from the dictionary
    ax.set_xlim(ExtentDict['Extent'][0],ExtentDict['Extent'][1])
    ax.set_ylim(ExtentDict['Extent'][2],ExtentDict['Extent'][3])
    #   Label the axes.
    ax.set_xlabel(r"$\Delta$ RA ('')")
    ax.set_ylabel(r"$\Delta$ DEC ('')")
    #   Return the extent dictionary as it used for drawing a circle in the Moment0 panel.
    return ExtentDict
    
def AddCenterToMomMap(ax):
    """
        This function adds a black X mark to the center point of a moment map.  Since the maps are formatted in delta RA and delta DEC from a center point, this is always 0,0
    """
    ax.plot(0,0,marker='x',color='k',markersize=10)
    return ax
    
def AddArrowToMomMap(ax,Angle,Length):
    """
        This function adds an arrow to a moment map panel to represent the position angle
    """
    #   First set the center coordinates of the arrow as the center of the galaxy
    CentX=0.
    CentY=0.
    #   Make an empty pair of sizes
    dX_Opts=np.zeros(2)
    dY_Opts=np.zeros(2)
    #   Get the distances to the edges of the plot from the center in both X and Y
    for i in range(2):
        dX_Opts[i]=np.abs(CentX-ax.get_xlim()[i])
        dY_Opts[i]=np.abs(CentY-ax.get_ylim()[i])
    #   Set the length of the arrow to be the minimum of all the distances calculated in order to ensure it ends inside the panel, but extends for a reasonable length
    L=1.*np.min((np.min(dX_Opts),np.min(dY_Opts)))
    #   Use the length and the angle to get the length in x and y for the arrow
    dX=-L*np.cos((Angle+90.)*np.pi/180.)
    dY=L*np.sin((Angle+90.)*np.pi/180.)
    #   Set the arrow width
    awidth=1.
    #   Add the arrow to the panel
    ax.arrow(CentX,CentY,dX,dY,ls=':',color='k',shape='full',width=awidth,head_width=10*awidth)
    return ax
    
def MakeMomMap(Cube,Moment,MaskSwitch):
    """
        This function constructs a moment 0, moment 1, or moment 2 array from a cube.
    """
    #       Get the map size
    Size=(np.shape(Cube['Data'])[1],np.shape(Cube['Data'])[2])
    #       Initialize the moment map and flux tot arrays
    MomMap=np.zeros(Size)
    FTot=np.zeros(Size)
    #   Decide whether to use the full data cube or the masked data cube.
    if MaskSwitch == 0:
        DUse=Cube['Data']
    else:
        DUse=Cube['MaskedData']
    #   The moment 0, 1, and 2 maps all need to start with a total flux map
    for i in range(np.shape(Cube['Data'])[1]):
        for j in range(np.shape(Cube['Data'])[2]):
            FTot[i,j]=np.nansum(DUse[:,i,j])
    #       When using unmasked data, the low flux cells can cause problems, so artificially set them to zero
    for i in range(np.shape(Cube['Data'])[1]):
        for j in range(np.shape(Cube['Data'])[2]):
            if FTot[i,j] <1.e-7:
                FTot[i,j]=0.
    #   Now construct the moment maps
    if Moment == 0:
        #   The moment 0 map is simply the total flux map.
        MomMap=FTot
    elif Moment==1 or Moment==2:
        #   For the moment 1 or moment 2 maps, we need to get the flux weighted average velocity in each pixel
        #       Loop through all the pixels
        for i in range(np.shape(Cube['Data'])[1]):
            for j in range(np.shape(Cube['Data'])[2]):
                #Loop through the channels
                for k in range(np.shape(Cube['Data'])[0]):
                    #   Get the flux-weighted sum for each channel
                    MomMap[i,j]+=(Cube['CubeVels'][k]/1000.)*DUse[k,i,j]
            #   Set the pixels with no flux to NaN's
            if FTot[i,j] ==0.:
                FTot[i,j]=FTot[i,j]/0.
        #   Normalize the map by the fluxes to get the moment 1 map
        MomMap=MomMap/FTot
    #   For the moment 2 map, go back through the pixels one more time
    if Moment ==2:
        #   Copy the moment 1 map to a temporary array
        VAvgMap=copy.copy(MomMap)
        #   Reset the moment map array to zeros
        MomMap=np.zeros(Size)
        #   Loop through all pixels and channels
        for i in range(np.shape(Cube['Data'])[1]):
            for j in range(np.shape(Cube['Data'])[2]):
                for k in range(np.shape(Cube['Data'])[0]):
                    #   Sum up the flux weighted squared difference
                    MomMap[i,j]+=(Cube['CubeVels'][k]/1000.-VAvgMap[i,j])**Moment*DUse[k,i,j]
        #   Now get the moment 2 map by taking the root of the mean square differences
        MomMap=np.sqrt(MomMap/FTot)
    #   Return either the moment 0, moment 1, or moment 2 map
    return MomMap
    
def AddVelContoursToMomentPlot(ax,Cube,VSys,ObsMap,XC,YC,Model):
    """
        This function adds velocity contours to a moment 1 map panel from some other cube
    """
    #   First make a moment map for the cube that will be used
    MomMap=MakeMomMap(Cube,1,0)
    #   Next get the extent of the panel using the center and the cube that the contours will be from
    ExtentDict=DetermineExtent(XC,YC,Cube)
    
    #   Figure out the velocity width of the cube
    #       Get the outermost velocity
    VOut=Model['VROT'][-1]
    #       Use the model to get the inclination corrected outermoste velocity
    VSinI=VOut*np.sin(Model['INCLINATION'][0]*np.pi/180.)
    #       Set the model profile width
    VWidth=2.*VSinI
    #   Set dV to be twice this width
    dV=1.0*VWidth
    #       Set the contour line type and linewidth
    lTypes=('-')
    LW=1
    #   Set an array of a X and Y corresponding to the pixels using the extent dictionary
    X=np.linspace(ExtentDict['Extent'][2],ExtentDict['Extent'][3],ExtentDict['nY'])
    Y=np.linspace(ExtentDict['Extent'][0],ExtentDict['Extent'][1],ExtentDict['nX'])
    #   Turn these into a meshgrid for use with the contour function
    YY2,XX2=np.meshgrid(Y,X)
    #   Set the velocity steps for the contours
    delV=dV/5.
    #   Set the contour levels
    CLevels=np.zeros(11)
    j=0
    for i in range(-5,6):
        CLevels[j]=VSys+i*dV/7.
        j+=1
    #   Draw on the contours
    ax.contour(YY2,XX2,MomMap,colors='magenta',linewidths=LW,linestyles=lTypes,levels=CLevels)
    #       Add in a thicker contour with vsys
    CLevels=np.array([VSys])
    ax.contour(YY2,XX2,MomMap,colors='magenta',linewidths=LW*2.5,linestyles=lTypes,levels=CLevels)
    #   Return the panel
    return ax
    
