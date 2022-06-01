import numpy as np
from decimal import Decimal

import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator

from matplotlib.patches import Ellipse

"""
    This module contains routines to make the full fit comparison diagnostic plot for a given galaxy.  It contains the routines:
    MakeComparisonPlot --> This function makes the full comparison plot.
    NameComparisonPlot --> This function names the particular diagnostic plot.
    AddObjectLabels --> This function adds a number of useful labels to the plot that give some basic measurements.
    AddGoodnessofFitLabels --> This function adds a set basic goodness fit measurements to the plot.
"""


def MakeComparisonPlot(ObjDict,FolderDict,CubeInfo,FitParams,FittingOptions,AnalysisFncs):
    """
        This function makes a diagnostic plot that compares a set of different Tilted Ring dictionaries (placed in an array).
    """
    #   Start by setting the plot parameters to some useful defaults
    BasePlotParams={'font.size': 15,'axes.linewidth':2
        ,'xtick.major.size':6,'xtick.minor.size':3
        ,'xtick.major.width':1,'xtick.minor.width':1
            ,'ytick.major.size':6,'ytick.minor.size':3
            ,'ytick.major.pad':10,'xtick.major.pad':10
            ,'ytick.major.width':1,'ytick.minor.width':1
            ,'xtick.labelsize':15 ,'ytick.labelsize':15
            ,'axes.labelsize': 20
            ,'legend.fontsize': 18
                }
    #   Update the default parameters to the new ones.
    matplotlib.rcParams.update(BasePlotParams)
    #   Name the plot
    PltName=NameComparisonPlot(ObjDict,FolderDict)
    #   Open the plot
    Fw=10.
    Fh=10.
    fig=plt.figure(figsize=(Fw,Fh))
    #   Make all the quick measurements from the cube -- they will be used for labels and in some subplots
    CubeMeasureDict=AnalysisFncs['CubeFnc']['MakeQuickCubeMeasures'](ObjDict,CubeInfo,AnalysisFncs)
    #   Write on all the labels due to the cube fit
    yTextLoc=AddObjectLabels(ObjDict,fig,CubeMeasureDict)
    #   Add the goodness of fit labels
    AddGoodnessofFitLabels(fig,FitParams,yTextLoc,FittingOptions,AnalysisFncs['PlotFncs'])
    #   Set the placement parameters that are used to organize all the different panels in the diagnostic plot
    base=-0.1
    left=0.1
    w=0.45
    h=w/2
    buf=0.15
    #   Make a major axis PV Plot
    #       Set the location of the plot
    placement=[left-0.8*(w+buf),base+2*(h+buf),w*0.8,h]
    #       Make the greyscale PV plot (see DiagnosticPlots/PVPlotFncs.py)
    PV_ax=AnalysisFncs['PlotFncs']['BasePVPlot'](fig,placement,CubeInfo,np.abs(CubeInfo['CubeHeader']['CDELT1'])*3600.)
    #       Add the RC's for each fit to the PV plot
    AnalysisFncs['PlotFncs']['AddRCs_to_PVPlot'](PV_ax,FitParams,FittingOptions,CubeInfo)
    #   Now add the moment map plots
    #       Start by placing the center point as the cube header's reference pixel (this will be set to zero in the moment map panels).
    XC=CubeInfo['CubeHeader']['CRPIX1']
    YC=CubeInfo['CubeHeader']['CRPIX2']
    #       Set the location of the moment 0 map
    placement=[left-0.8*(w+buf),base+1*(h+buf),w,h]
    #       Make the moment 0 map (see DiagnosticPlots/MomentMapPlotFncs.py)
    Moment=0
    AnalysisFncs['PlotFncs']['MomentPlot'](fig,placement,Moment,ObjDict,CubeInfo,XC,YC)
    #       Set the location of the moment 1 map
    placement=[left-0.8*(w+buf),base+0.*(h+buf),w,h]
    #       Make the moment 1 map
    Moment=1
    AnalysisFncs['PlotFncs']['MomentPlot'](fig,placement,Moment,ObjDict,CubeInfo,XC,YC)
    #       Set the location of the moment 2 map
    placement=[left-0.8*(w+buf),base-1.*(h+buf),w,h]
    #       Make the moment 2 map
    Moment=2
    AnalysisFncs['PlotFncs']['MomentPlot'](fig,placement,Moment,ObjDict,CubeInfo,XC,YC)
    #   Add the keyword plots (see DiagnosticPlots/DiagnosticPlotFuncs.py)
    #       Do the rotation velocity plot
    placement=[left,base+2*(h+buf),w,h]
    Key="VROT"
    AnalysisFncs['PlotFncs']['KeywordPlot'](fig,placement,Key,ObjDict,FitParams,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
     #       Do the surface density plot
    placement=[left+1*(w+buf),base+2*(h+buf),w,h]
    Key="SURFDENS"
    AnalysisFncs['PlotFncs']['KeywordPlot'](fig,placement,Key,ObjDict,FitParams,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
    #       Do the inclination plot
    placement=[left,base+1*(h+buf),w,h]
    Key="INCLINATION"
    AnalysisFncs['PlotFncs']['KeywordPlot'](fig,placement,Key,ObjDict,FitParams,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
    #       Do the position angle plot
    placement=[left+(w+buf),base+1*(h+buf),w,h]
    Key="POSITIONANGLE"
    AnalysisFncs['PlotFncs']['KeywordPlot'](fig,placement,Key,ObjDict,FitParams,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
    #       Do the systemic velocityplot
    placement=[left,base-1*(h+buf),w,h]
    Key="VSYS"
    AnalysisFncs['PlotFncs']['KeywordPlot'](fig,placement,Key,ObjDict,FitParams,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
    #       Do the x center plot
    placement=[left,base-0*(h+buf),w,h]
    Key="XCENTER"
    AnalysisFncs['PlotFncs']['KeywordPlot'](fig,placement,Key,ObjDict,FitParams,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
    #       Do the y center plot
    placement=[left+(w+buf),base-0*(h+buf),w,h]
    Key="YCENTER"
    ax1=AnalysisFncs['PlotFncs']['KeywordPlot'](fig,placement,Key,ObjDict,FitParams,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
    #       Add a legend to the plot for easier reading
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels,bbox_to_anchor=(1.05, -0.5), loc='upper left', borderaxespad=0)
    #   Save the plot
    plt.savefig(PltName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()

def NameComparisonPlot(ObjDict,FolderDict):
    """
        Thsi function names the plot using the target fit folder and the WALLABY detection name.
    """
    PltName=FolderDict['FitsComparisionFolder']+ObjDict['NumStr']+"_FitComparisonPlot.png"
    return PltName


def AddObjectLabels(ObjDict,fig,CubeMeasureDict):
    """
        This function adds a number of labels to the plot
    """
    #   First add the name of the galaxy to the plot
    fig.text(.6,.95, ObjDict['CatEntry']['name'] , ha='center',rotation=0,va='center',size=27)
    #   Next add the size of the major axis from SoFiA in beams
    yTextLoc=0.90
    LabelStr="Ell_Maj \t=\t".expandtabs()+str(round(CubeMeasureDict['Diameter'],3))+" beams"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    #   Then add the 'SoFiA radius' by dividing Ell_Maj by 2
    yTextLoc-=0.1
    LabelStr="Ell_Maj/2 \t=\t".expandtabs()+str(round(CubeMeasureDict['Radius'],3))+" '' "
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    #   Add the predicted value of R_HI from the HI mass-size relation for comparison with the SoFiA measurements.
    yTextLoc-=0.1
    LabelStr="RHI \t\t=\t".expandtabs()+str(round(CubeMeasureDict['RHI'],3))+" '' "
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    #   Add the estimated mass (calculated from the total cube flux, the redshift velocity and the Hubble flow distance estimate).
    yTextLoc-=0.1
    MStr='%2E'%Decimal(CubeMeasureDict['Mass'])
    LabelStr="Mass\t\t= \t".expandtabs()+MStr+r" M$_{\odot}$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    #   Add the distance estimated using the redshift and the Hubble flow.
    yTextLoc-=0.1
    LabelStr="Distance \t=\t".expandtabs()+str(round(CubeMeasureDict['Distance'],3) )+" Mpc"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=25,color='k')
    yTextLoc-=0.15
    #   Return the current height of the labels
    return yTextLoc

def AddGoodnessofFitLabels(fig,FitParams,yTextLoc,FittingOptions,PlotFncs):
    """
        This function adds the basic goodness of fit statistics calculated with loading in and analyzing the various fits.
    """
    #   Add a Goodness of fit heading
    fig.text(1.45,yTextLoc, r"Goodness of Fit" , ha='center',rotation=0,va='center',size=22)
    yTextLoc-=0.1
    #   Now loop through all fits
    for i in range(FittingOptions['nTotFits']):
        #   Get the line color from the plot functions
        linecol=PlotFncs['ColorSelection'](i)
        #   Get the line label
        linelabel=FitParams[i]['Label']
        #   If there is a successful fit, add in the reduced Chi^2 value to the plot
        if FitParams[i]['FITAchieved']:
            chi2=FitParams[i]['CHI2']
            LabelStr=linelabel+"\t = \t".expandtabs()+str(round(chi2,3) )
            fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=18,color=linecol)
            yTextLoc-=0.05


    
    

