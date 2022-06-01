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

from . import MomentMapPlotFncs as MMPF
from . import PVPlotFncs as PVPF
from . import GeneralPlotDefs as GPD

"""
    This module contains a number of plotting functions that are useful for both diagnostic and average model plots.  It contains the routines:
    PlotFnc --> This contains a dictionary of functions that can be used by other modules without directly importing this module.
    KeywordPlot --> This function makes a panel showing the parameter values for a variety of different models for a specific parameter (RC, SD, center, orientation, etc.).
    KeywordPlot_SingleModel --> This function makes a panel like the KeywordPlot routine but only deals with a single model.
    
    SetXY --> This selects the X and Y values for a keyword panel line
    AddErrorBars --> For certain keyword panels, this adds the uncertainties to the plotted lines.
    LabelSelection --> This function sets the labels for the keyword plots
    GetYLim --> This function sets the upper limits on the keyword plots
    KeywordPlotFmt --> This function formats a keyword panel plot
"""



def PlotFnc():
    """
        This function stores other routines in a dictionary so that other modules can use them without directly importing this module.  It combines routines found in other modules for cleaner organization/calling
    """
    PltFncDict={'ColorSelection':GPD.ColorSelection,'MomentPlot':MMPF.MomentPlot,'AddCenterToMomMap':MMPF.AddCenterToMomMap,'AddArrowToMomMap':MMPF.AddArrowToMomMap,'AddVelContoursToMomentPlot':MMPF.AddVelContoursToMomentPlot,'MakeMomMap':MMPF.MakeMomMap,'KeywordPlot':KeywordPlot,'BasePVPlot':PVPF.BasePVPlot,'AddRCs_to_PVPlot':PVPF.AddRCs_to_PVPlot,'AddPVContoursToPlot':PVPF.AddPVContoursToPlot,'AddSingleRC_to_PVPlot':PVPF.AddSingleRC_to_PVPlot,'AddCentLinesToPVPlot':PVPF.AddCentLinesToPVPlot,'AvgModelPVPlot':PVPF.AvgModelPVPlot,'KeywordPlot_SingleModel':KeywordPlot_SingleModel}
    return PltFncDict
        
def KeywordPlot(fig,placement,Key,ObjDict,FitParams,CubeMeasureDict,FittingOptions,AstroFncs):
    """
        This function makes a plot showing all different model values for a given parameter indicated by the keyword
    """
    #   Set the type of diagnostic plot this is
    DiagnosticSwitch=0
    #   Draw the axis
    ax=fig.add_axes(placement)
    #   Set line and marker sizes
    LW=1
    MW=10
    #   Set the fit counters to zero
    Step=0
    Advance=0
    #   Determine if a catalogue line should be drawn
    drawCatalogueLine,CatVal=CatLine(Key,ObjDict,AstroFncs,CubeMeasureDict)
    if drawCatalogueLine:
        ax.axhline(y=CatVal,color='k',ls='--',lw=LW)
    #   Set Ymax to zero for setting axis limits
    YMax=0.
    #   Now loop through all fits
    for i in range(FittingOptions['nTotFits']):
        #   Get the line color
        linecol=GPD.ColorSelection(i)
        #   Get the line label
        linelabel=FitParams[i]['Label']
        #   If there is a successful fit, start the plotting steps
        if FitParams[i]['FITAchieved']:
            #   Set the X and Y values
            X,Y=SetXY(Key,FitParams[i])
            #   Adjust the max value of Y
            if np.max(Y) > YMax:
                YMax=np.max(Y)
            #   Plot the X-Y points
            ax.plot(X,Y,marker='.',ls=':',color=linecol,lw=LW,markersize=MW,label=linelabel)
            #   For some key words, add in error bars
            AddErrorBars(ax,LW,MW,linelabel,linecol,Key,X,Y,FitParams[i],DiagnosticSwitch)
    #   After the loop, add vertical lines for the radius
    ax.axvline(x=CubeMeasureDict['Radius'], color='k',ls='--',lw=LW)
    ax.axvline(x=CubeMeasureDict['RHI'], color='k',ls='-.',lw=LW)
    #   Finally format the plot
    KeywordPlotFmt(ax,Key,YMax,CatVal,DiagnosticSwitch)
    return ax

def KeywordPlot_SingleModel(fig,placement,Key,ObjDict,Model,CubeMeasureDict,FittingOptions,AstroFncs):
    """
        This function makes a plot showing a specific model's values for a given parameter indicated by the keyword
    """
    #   Set the type of diagnostic plot this is
    DiagnosticSwitch=1
    #   Draw the axis
    ax=fig.add_axes(placement)
    #   Set line and marker sizes
    LW=1
    MW=10
    #   Set the fit counters to zero
    Step=0
    Advance=0
    #   Determine if a catalogue line should be drawn
    drawCatalogueLine,CatVal=CatLine(Key,ObjDict,AstroFncs,CubeMeasureDict)
    if drawCatalogueLine:
        ax.axhline(y=CatVal,color='k',ls='--',lw=LW)
    #   Set Ymax to zero for setting axis limits
    YMax=0.
    #   Now loop through all fits
        #   Get the line color -- this sets it to black
    linecol=GPD.ColorSelection(-1)
        #   Get the line label -- for a single model, we don't want a generally want label
    linelabel=" "
        #   Set the X and Y values
    X,Y=SetXY(Key,Model)
        #   Adjust the max value of Y
    if np.max(Y) > YMax:
        YMax=np.max(Y)
        #   Plot the X-Y points
    LLabel="Measured"
    L1,=ax.plot(X,Y,marker='.',ls='-',color=linecol,lw=LW,markersize=MW,label=LLabel)
        #   For some key words, add in error bars
    AddErrorBars(ax,LW,MW,LLabel,linecol,Key,X,Y,Model,DiagnosticSwitch)
    #   For the surface density, we want to also show the face-on values
    if Key == 'SURFDENS':
        #   Set Y to the face-on SD
        YSpec=Model['SURFDENS_FACEON']
        LLabel="Face On"
        L2,=ax.plot(X,YSpec,marker='.',ls=':',color='grey',lw=LW,markersize=MW,label=LLabel)
        #   Add a lagend to the SD panel
        ax.legend(handles=[L1,L2],prop={'size': 15})

    KeywordPlotFmt(ax,Key,YMax,CatVal,DiagnosticSwitch)
    return ax

def SetXY(Key,FitParams):
    """
        This function sets the X and Y values for plotting
    """
    #   Start by seting to the bass X and Y
    X=FitParams['R']
    Y=FitParams[Key]
    #   If the key is the surface density, adjust to using the 'R_SD' for X
    if Key == 'SURFDENS':
        X=FitParams['R_SD']
    return X,Y
    
def AddErrorBars(ax,LW,MW,linelabel,linecol,Key,X,Y,FitParams,DiagnosticSwitch):
    """
        This function adds error bars to the keyword plots if appropriate.  Different plot styles use different sets of error bars
    """
    #   Start by assuming no error bars
    PlotErr=False
    #   For model comparison diagnostic plots, show the error bars for Vrot, Inc, and PA
    if DiagnosticSwitch ==0 :
        if Key == "VROT" or Key == "INCLINATION" or Key == "POSITIONANGLE":
            PlotErr=True
    #   For average model diagnostic plots, show the error bars for Vrot and SD
    elif DiagnosticSwitch == 1:
        if Key == "VROT" or Key == "SURFDENS" or Key == "SURFDENS_INC":
            PlotErr=True
    if PlotErr :
        #   If adding errorbars, set the keyword error
        KeyErr=Key+"_ERR"
        #   Get the errorbars and set them absolutes
        YErr=FitParams[KeyErr]
        YErr=np.abs(YErr)
        #   Add in the error bars
        ax.errorbar(X,Y,yerr=YErr,marker='.',ls='-',color=linecol,lw=LW,markersize=MW,label=linelabel)

def LabelSelection(Key,DiagnosticSwitch):
    """
        This function sets the y-labels for the various keyword panels in the diagnostic plots.
    """
    if Key == "VROT":
        #   This can be used for either normal diagnostic or single model plots with different labels for VROT and SURFDENS
        if DiagnosticSwitch == 0:
            ylabel=r"$V_{\phi}$"
        elif DiagnosticSwitch== 1:
            ylabel =r"VROT_kin (km/s)"
    elif Key == "SURFDENS":
        if DiagnosticSwitch == 0:
            ylabel=r"$\Sigma$ "
        elif DiagnosticSwitch== 1:
            ylabel = r"SD_kin (M$_{\odot}$ pc$^{-2}$)"
    elif Key == "INCLINATION":
        ylabel=r"$I$ "
    elif Key == "POSITIONANGLE":
        ylabel=r"$PA$ "
    elif Key == "VSYS":
        ylabel=r"$V_{sys}$ "
    elif Key == "VDISPERSION":
        ylabel=r"$\sigma_{r}$ "
    elif Key == "XCENTER":
        ylabel=r"$X$ "
    elif Key == "YCENTER":
        ylabel=r"$Y$ "
    else:
        ylabel=" "
    return ylabel

def GetYLim(Ymax,LineVal,Key,DiagnosticSwitch):
    """
        This function sets the upper limit for a keyword plot
    """
    #   This can be used for either normal diagnostic or model average plots with different surface density units
    #   Generally set the upper limit to the maximum value in set of models
    YLim=Ymax
    #   For Vrot, make sure there's a bit of upper room, and set a hard limit at 400 km/s
    if Key == "VROT":
        YLim=np.min((1.1*YLim,400.))
    #   For SD, the units are different in the fit comparison and model average plots, so set the limits appropriately.
    elif Key == "SURFDENS":
        if DiagnosticSwitch== 0:
            YLim=np.min((1.1*YLim,1.e-3))
        elif DiagnosticSwitch== 1:
            YLim=np.min((1.1*YLim,5.e1))
    elif Key == "INCLINATION":
        YLim=np.max((YLim+2.,LineVal+2))
    elif Key == "POSITIONANGLE":
        YLim=np.max((YLim+2.,LineVal+1))
    elif Key == "VSYS":
        YLim=np.max((YLim,LineVal))+2
    elif Key == "VDISPERSION":
        YLim=np.min((YLim,30.))
    elif Key == "XCENTER":
        YLim=np.max((YLim,LineVal))+0.5
    elif Key == "YCENTER":
        YLim=np.max((YLim,LineVal))+0.5
    else:
        YLim=np.min(YLim)
    return YLim

def CatLine(Key,ObjDict,AstroFncs,CubeMeasureDict):
    """
        This function draws a horizontal line from a catalogue value for certain keyword plots.  This is useful for seeing how the kinematic measurements compare to the catalogue values.
    """
    #   Assume that a catalogue value won't be drawn and make an empty variable to return
    drawCatVal=False
    CatVal=None
    #   For inclination, PA, Vsys, X, and Y, get the catalogue measurements in the appropriate units.
    if Key == "INCLINATION":
        drawCatVal=True
        #   The SoFiA catalogue ell_min and ell_maj values need to be converted to an inclination via cos(i)=ell_min/ell_maj
        CatVal=np.arccos(ObjDict['CatEntry']['ell_min']/ObjDict['CatEntry']['ell_maj'])*180./np.pi
    elif Key == "POSITIONANGLE":
        drawCatVal=True
        #   For PA, the SoFiA values differ from the kinematics by 180 degrees
        CatVal=ObjDict['CatEntry']['kin_pa']+180.
        #   Make sure the value is between 0<PA<360
        if CatVal > 360.:
            CatVal=CatVal-360.
    elif Key == "VSYS":
        #   For Vsys, the catalogue frequency value needs to be converted to a velocity
        drawCatVal=True
        Freq=ObjDict['CatEntry']['freq']
        RestFreq=1.42040575179E+09
        CatVal=AstroFncs['RedShiftConv'](Freq,RestFreq)/1000.
    elif Key == "XCENTER":
        drawCatVal=True
        CatVal=CubeMeasureDict['CentPix'][0]
    elif Key == "YCENTER":
        drawCatVal=True
        CatVal=CubeMeasureDict['CentPix'][1]
    return drawCatVal,CatVal

def KeywordPlotFmt(ax,Key,YMax,CatVal,DiagnosticSwitch):
    """
        This function formats a diagnostic plot nicely.
    """
    #   Label the plot
    ax.set_xlabel(r"R ('')")
    ylabel=LabelSelection(Key,DiagnosticSwitch)
    ax.set_ylabel(ylabel)
    #   Set the Ymax
    YMax=GetYLim(YMax,CatVal,Key,DiagnosticSwitch)
    ax.set_ylim(top=YMax)
    #   Set the lower limit to 0 for Vrot and SD
    if Key == "VROT" or Key == "SURFDENS":
        ax.set_ylim(bottom=0.)
    #   Add minor tick marks to the panel.
    ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n=5))
