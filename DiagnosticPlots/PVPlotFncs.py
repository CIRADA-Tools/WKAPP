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

from . import GeneralPlotDefs as GPD

"""
    This module contains a number of plotting functions that are useful for both diagnostic and average model plots.  IThis group is focus on position-velocity plot functions. It contains the routines:

    BasePVPlot --> This function makes a basic position-velocity plot panel for some data.
    AddRCs_to_PVPlot --> This function adds a set of inclination corrected rotation curves to a PV plot panel.
    AddSingleRC_to_PVPlot --> This function adds a single RC to a PV panel.
    DrawSingleRC_onPVPanel --> This function draws the RC from a particular model on a panel.  It is used by both AddRCs_to_PVPlot and AddSingleRC_to_PVPlot.
    AddPVContoursToPlot --> This function adds contours to a PV plot based on some other cube.
    AvgModelPVPlot --> This makes a PV diagnostic panel for the average model diagnostic plots
    GetPVNoise --> This function gets the noise in some PV diagram array
    SetCornerLims --> This function gets the limits to sum over that represent the 4 corners of the PV diagram
    CornerSum --> This function sums up the quantities in a corner of an array.
"""

def BasePVPlot(fig,placement,CubeInfo,PixSize):
    """
        This function makes a PV panel plot
    """
    #   First draw the panel onto the larger canvas
    ax=fig.add_axes(placement)
    #   First set an array of X in 'pixels' corresponding to the already made PV array
    X=np.linspace(1,np.shape(CubeInfo['PV'])[0],np.shape(CubeInfo['PV'])[0])
    #   Set the mid point to 0
    XMid=np.shape(CubeInfo['PV'])[0]/2
    X=X-XMid
    #   Adjust the position coordinates into arcseconds
    X=X*PixSize
    #   Set the velocity array to be the channel velocities.
    V=CubeInfo['CubeVels']/1000.
    #   Make a meshgrid
    VV,XX=np.meshgrid(V,X)
    #   Use the cube PV array to make a colormesh PV panel that is grey-scaled
    ax.pcolormesh(XX,VV,CubeInfo['PV'],cmap='Greys',shading='auto')
    #   Label the Y-axis -- don't label the X-axis as it may be a major or minor axis PV diagram
    ax.set_ylabel(r"V (km/s)" )
    #   Format the panel with minor tick marks
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=4))
    #   Return the full panel
    return ax
    
def AddCentLinesToPVPlot(ax,VSys,XCenter):
    """
        This function adds center lines to a PV panel through the midpoint and systemic velocity
    """
    ax.axhline(y=VSys,ls=':',color='green',linewidth=2)
    ax.axvline(x=XCenter,ls=':',color='green',linewidth=2)
    return ax
    
def AddPVContoursToPlot(ax,CubeInfo,PixSize,PVNoise):
    """
        This function adds contours to an existing PV diagram from a different cube.  The contours levels are based on the noise in the existing PV panel.
    """
    #   Set the contour levels at 1,3, and 5 times the noise
    CLevels=np.array([1.,3,5.])*PVNoise

    #   Set different line types for each contour
    lTypes=(':','--','-')
    #   Set the thickness of the contours
    LW=1
    #   Set a velocity array from the channel velocities
    V2=CubeInfo['CubeVels']/1000.
    #   Set an array of X based on the size of the cube PV diagram
    X=np.linspace(1,np.shape(CubeInfo['PV'])[0],np.shape(CubeInfo['PV'])[0])
    #   Reset the X array to place the midpoint at 0
    XMid=np.shape(CubeInfo['PV'])[0]/2
    X=X-XMid
    #   Convert from pixels to arcseconds
    X=X*PixSize
    #   Construct a meshgrid needed for the contours
    VV2,XX2=np.meshgrid(V2,X)
    #   Draw the contours onto the plot in magenta.
    ax.contour(XX2,VV2,CubeInfo['PV'],levels=CLevels,colors='magenta',linewidths=LW,linestyles=lTypes)
    #   Return the panel
    return ax
    
def AddRCs_to_PVPlot(ax,FitParams,FittingOptions,CubeInfo):
    """
        This function draws a set of rotation curves onto a PV plot.  This is used in the model comparison diagnsotic plot.
    """
    #   Set the line and marker sizes
    LW=1
    MW=10
    #   Set the pixel size to get the points
    pixSize=np.abs(CubeInfo['CubeHeader']['CDELT2'])*3600.
    #       Loop through all the models
    for i in range(FittingOptions['nTotFits']):
        #   Get the line color
        linecol=GPD.ColorSelection(i)
        #   Get the line label
        linelabel=FitParams[i]['Label']
        #   If there is a successful fit, start the plotting steps
        if FitParams[i]['FITAchieved']:
            #   Draw the RC projected onto the PV diagram for each model
            DrawSingleRC_onPVPanel(ax,FitParams[i],linecol,LW,MW,linelabel)

            
def AddSingleRC_to_PVPlot(ax,Model,CubeInfo):
    """
        This function adds a single rotation curve to a PV panel for the average model diagnostic plots
    """
    #   Set the linewidth, marker size, color, and a general label
    LW=1
    MW=10
    linecol='blue'
    linelabel=" "
    #   Draw the RC from the model onto the panel
    DrawSingleRC_onPVPanel(ax,Model,linecol,LW,MW,linelabel)
 
def DrawSingleRC_onPVPanel(ax,Model,linecol,LW,MW,linelabel):
    """
        This function draws the RC from a specific model onto a given PV panel.
    """
    #   Do the positive radius portion
    #       Set the X values
    X1=+   Model['R']
    #       Set the Y values using Vsys + Vrot*sin(Inc)
    V1=Model['VSYS']+(Model['VROT']*np.sin(Model['INCLINATION']*np.pi/180.))
    #   Plot the positive portion
    ax.plot(X1,V1,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)
    #   Do the same thing for the negative radius portion
    #       Set the X values of the line
    X2=-Model['R']
    #       Set the Y values using Vsys - Vrot*sin(Inc)
    V2=Model['VSYS']-(Model['VROT']*np.sin(Model['INCLINATION']*np.pi/180.))
    #   Plot the negative portion
    ax.plot(X2,V2,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)





def AvgModelPVPlot(fig,placement,ObjDict,AvgModel,CubeInfo,FittingOptions,AnalysisFncs,MajorMinorSwitch,ModelCube):

    """
        This function makes the PV plot for a single 'average' model
    """
    #   First use the average model to find the center to be used for calculating the PV diagram
    CentPix=[AvgModel['XCENTER'][0],AvgModel['YCENTER'][0]]
    #   Now set the Position angle and center for the PV plot
    if MajorMinorSwitch==0:
        PAUse=AvgModel['POSITIONANGLE'][0]
        XLabel=r"Major Axis ('')"
        Center=AvgModel['XCENTER'][0]
    elif MajorMinorSwitch ==1:
        PAUse=AvgModel['POSITIONANGLE'][0]+90.
        XLabel=r"Minor Axis ('')"
        Center=AvgModel['YCENTER'][0]
    #   Normalize the useable position angle
    if PAUse > 360.:
        PAUse=PAUse-360.
    elif PAUse < 0.:
        PAUse=PAUse+360.
    #   The PV routine needs the beamsize in pixels to figure out how things should be cut
    BeamSize_Pix=CubeInfo['CubeHeader']['BMAJ']/np.abs(CubeInfo['CubeHeader']['CDELT1'])
    #   Construct the PV diagram for the data cube using the average geometry
    PV=AnalysisFncs['CubeFnc']['ConstructModelBasedPVDiagram'](CubeInfo['Data'],PAUse,AvgModel,BeamSize_Pix,CubeInfo)
    #   Add this all to a PV dictionary
    PVDict={'PV':PV,'CubeVels':CubeInfo['CubeVels']}
    #   Get the noise for the PV diagram
    PVNoise=GetPVNoise(PV)
    #   Get the pixel size in arcseconds
    PixSize=np.abs(CubeInfo['CubeHeader']['CDELT1']*3600.)
    #   Make the grey-scale PV plot
    ax=BasePVPlot(fig,placement,PVDict,PixSize)
    #       Add lines for vsys and the center point to the PV diagram
    ax=AddCentLinesToPVPlot(ax,AvgModel['VSYS'][0],0.)
    #   Make a PV diagram of the model cube
    ModelPV=AnalysisFncs['CubeFnc']['ConstructModelBasedPVDiagram'](ModelCube['Data'],PAUse,AvgModel,BeamSize_Pix,CubeInfo)
    #   Adjust the PV dictionary to use the new PV array
    PVDict={'PV':ModelPV,'CubeVels':ModelCube['CubeVels']}
    #   Add the model PV diagram contours overtop the observed PV map
    ax=AddPVContoursToPlot(ax,PVDict,PixSize,PVNoise)
    #   Set the xlabel for the plot
    ax.set_xlabel(XLabel)
    return ax


def GetPVNoise(PV):
    """
        This function gets the noise for a PV plot using the 4 corners of the plot.  This is used to set the contour levels for another model overplotted on the main plot
    """
    #   Set the relative size of the corners -- this must be large enough to be more than just zeros (which can occur due to the equal radial size extending beyond the pixel edges on a particular side)
    CornerSize=0.25
    #   Get the shape of the PV diagram
    Shape=np.shape(PV)
    #   Set a counter for the number of pixels and set the RMS to zero
    nPix=0
    RmsTot=0.
    #   Loop through each corner, setting the limits, and getting the squared flux sum and number of pixels
    for i in range(2):
        #   The i dimension is the X axis
        dimSwitch=0
        #   Set the limits on X for the particular corner
        XCornerLims=SetCornerLims(CornerSize,Shape,i,dimSwitch)
        for j in range(2):
            dimSwitch=1
            #   Set the limits on Y for the particular corner
            YCornerLims=SetCornerLims(CornerSize,Shape,j,dimSwitch)
            #   Make an array giving the X and Y limits for the corner together
            FullLims=[XCornerLims,YCornerLims]
            #   Sum up the number of pixels and the square of the flux in this corner
            nPix,RmsTot=CornerSum(PV,FullLims,nPix,RmsTot)
    #   Set the noise as the RMS of pixels in the 4 corners.
    PVNoise=np.sqrt(RmsTot/nPix)
    #   Return the noise
    return PVNoise

def SetCornerLims(CornerSize,Shape,step,dimSwitch):
    """
        This function sets the limits in X or Y in units of pixels for a given corner.
    """
    #   Either do the left/bottom limits
    if step == 0:
        Low=0
        High=int(CornerSize*Shape[dimSwitch])
        #   Or do the top/right limits
    else:
        High=int(Shape[dimSwitch]-1)
        Low=int(High-CornerSize*Shape[dimSwitch])
    Lims=[Low,High]
    return Lims

def CornerSum(PV,FullLims,nPix,RmsTot):
    """
        This function gets the number of pixels and the square of the flux in a corner of an array
    """
    #   Loop through all X values
    for i in range(FullLims[0][0],FullLims[0][1]):
        #   Loop through all Y values
        for j in range(FullLims[1][0],FullLims[1][1]):
            #   If there's any flux, add in the square and increase the counter.  When it's zero, it's because it's at an  edge, so don't count it.
            if PV[i][j] !=0.:
                nPix+=1
                RmsTot+=PV[i][j]**2.
    return nPix,RmsTot
