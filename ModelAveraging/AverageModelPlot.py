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
from astropy.coordinates import Angle

"""
    This module contains the routines needed to make the average model diagnostic plot for a specific galaxy.  It contains the routines:
    MakeAvgModelPlot --> This makes and saves the diagnostic plot for the 'average' model.
    NameComparisonPlot --> This function names the plot
    AddObjectLabels --> This function adds a set of labels to the plot.
    AvgPVPlots --> This function makes the pair of PV diagrams needed for the comparison plot.
"""

def MakeAvgModelPlot(ObjDict,FolderDict,CubeInfo,AvgModel,FittingOptions,AvgOutputDict,AnalysisFncs):
    """
        This function makes the diagnostic plot for the average Tilted ring model.
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
    #   Load in the full resolution MCG model cube realization for the moment map and PV plot comparisons.
    ModelCube=AnalysisFncs['CubeFnc']['BasicCubeAnalysis'](ObjDict['FullResAvgModelCubeName'])
    #   Name the plot
    PltName=NameComparisonPlot(ObjDict,FittingOptions,AvgOutputDict)
    #   Open the plot
    Fw=10.
    Fh=10.
    fig=plt.figure(figsize=(Fw,Fh))
    #   Make all the quick measurements from the cube -- they will be used for labels and in some subplots
    CubeMeasureDict=AnalysisFncs['CubeFnc']['MakeQuickCubeMeasures'](ObjDict,CubeInfo,AnalysisFncs)
    #   Write on all the labels due to the cube fit
    yTextLoc=AddObjectLabels(ObjDict,fig,CubeMeasureDict,AvgModel)
    
    #   Set the placement parameters that are used to organize all the different panels in the diagnostic plot
    base=-0.1
    left=0.1
    w=0.45
    h=w/2
    buf=0.15
    #   Add the keyword plots (see DiagnosticPlots/DiagnosticPlotFuncs.py)
    #       Add the rotation curve panel
    placement=[left,base+2*(h+buf),w,h]
    Key="VROT"
    AnalysisFncs['PlotFncs']['KeywordPlot_SingleModel'](fig,placement,Key,ObjDict,AvgModel,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
    #       Add the surface density panel
    placement=[left+1*(w+buf),base+2*(h+buf),w,h]
    Key="SURFDENS"
    AnalysisFncs['PlotFncs']['KeywordPlot_SingleModel'](fig,placement,Key,ObjDict,AvgModel,CubeMeasureDict,FittingOptions,AnalysisFncs['AstroFncs'])
    #       Add the pair of PV diagram plots
    AvgPVPlots(fig,ObjDict,AvgModel,CubeInfo,FittingOptions,AnalysisFncs,left,base,w,h,buf)
    #   Rescale the sizes/boxes for the moment maps
    ScaleFactor=1.25
    wnew=w*ScaleFactor
    hnew=h*ScaleFactor
    left=left-(wnew-w)/2.
    base=base-(hnew-h)/1.5
    #   Make the moment 0 map
    placement=[left-0.*(w+buf),base+1.0*(h+buf),wnew,hnew]
    Moment=0
    #       For the moment maps, it's necessary to give it the model center in pixels to center the maps
    XC=AvgModel['XCENTER'][0]
    YC=AvgModel['YCENTER'][0]
    #   Call the moment map function to make the basic map
    ax,MomData=AnalysisFncs['PlotFncs']['MomentPlot'](fig,placement,Moment,ObjDict,CubeInfo,XC,YC)
    #   Add a marker to the center of the map
    ax=AnalysisFncs['PlotFncs']['AddCenterToMomMap'](ax)
    placement=[left+1.*(w+buf),base+1.*(h+buf),wnew,hnew]
    #   Make the moment 1 map
    Moment=1
    ax,MomData=AnalysisFncs['PlotFncs']['MomentPlot'](fig,placement,Moment,ObjDict,CubeInfo,XC,YC,AvgModel)
    #   Add the center to moment 1 map
    ax=AnalysisFncs['PlotFncs']['AddCenterToMomMap'](ax)
    #   Add the model velocity contours to the observed moment 1 map
    ax=AnalysisFncs['PlotFncs']['AddVelContoursToMomentPlot'](ax,ModelCube,AvgModel['VSYS'][0],MomData,XC,YC,AvgModel)
    #   Add an arrow to the PV diagram indicating the position angle of the model
    ax=AnalysisFncs['PlotFncs']['AddArrowToMomMap'](ax,AvgModel['POSITIONANGLE'][0],0.5*ObjDict['MeasureCat']['SIZE'].values[0]*30./5.)
    #   Save the plot
    plt.savefig(PltName, format='png',bbox_inches='tight')
    #   Close the plot
    plt.close()
    
def NameComparisonPlot(ObjDict,FittingOptions,AvgOutputDict):
    """
        Name the diagnostic plot file
    """
    PltName=AvgOutputDict['OutputDirName']+AvgOutputDict['OutputBaseName']+"_DiagnosticPlot.png"
    return PltName


def AddObjectLabels(ObjDict,fig,CubeMeasureDict,AvgModel):
    """
        This function adds a number of labels to the plot
    """
    #   First add the name of the plot
    TextSize=18
    fig.text(.6,.95, ObjDict['CatEntry']['name'] , ha='center',rotation=0,va='center',size=27)
    #   Set a top location
    yTextLoc=0.50
    #   Set size of the vertical steps for the labels
    YTextStep=0.06
    
    #   Convert the RA into hours:minutes:seconds
    RA_Test=Angle(AvgModel['RA'][0],u.deg)
    RA_HMS=RA_Test.hms
    #   Format the RA string
    RA_HStr=""
    for i in range(len(RA_HMS)):
        Val=RA_HMS[i]
        if i < 2:
            Val=round(Val,0)
            ValS=str(Val).split('.')[0]
        else:
            Val=round(Val,1)
            ValS=str(Val)
        RA_HStr+=ValS
        if i==0:
            RA_HStr+="h"
        elif i == 1:
            RA_HStr+="m"
        elif i ==2:
            RA_HStr+="s"
    #       Do the same for the uncertainty in RA and turn it into arcseconds
    RA_Err=Angle(AvgModel['RA_ERR'],u.deg)
    RA_Err_Str=RA_Err.to_string(unit=u.degree)
    Err_Str=RA_Err.dms
    Err_Str=str(round(Err_Str[2],1))
    #Err_Str.split('.')
    Err_Str=Err_Str+"''"
    #   Add the RA label to the plot
    LabelStr="RA_model \t\t=\t".expandtabs()+RA_HStr+" $\pm$ " +Err_Str
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep

    #   Convert the DEC from degrees to degrees:minutes:arcseconds format
    DEC_Test=Angle(AvgModel['DEC'][0],u.deg)
    DEC_dms=DEC_Test.signed_dms
    #   Format the DEC string
    DEC_DStr=""
    for i in range(len(DEC_dms)):
        Val=DEC_dms[i]
        if i ==0 :
            if Val < 0:
                ValS="-"
        elif i < 3:
            ValS=str(Val).split('.')[0]
        else:
            ValS=str(round(Val,1))
        DEC_DStr+=ValS
        if i ==1:
            DEC_DStr+=r"$^{\circ}$"
        elif i ==2:
            DEC_DStr+="'"
        elif i == 3:
            DEC_DStr+="''"
    #   Do the same thing for the DEC uncertaintites.
    DEC_Err=Angle(AvgModel['DEC_ERR'],u.deg)
    DEC_Err_Str=DEC_Err.to_string(unit=u.degree)
    Err_Str=DEC_Err.dms
    Err_Str=str(round(Err_Str[2],1))
    Err_Str=Err_Str+"''"
    #   Add the DEC string to the plot
    LabelStr="DEC_model \t=\t".expandtabs()+DEC_DStr+" $\pm$ " +Err_Str
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    #   Add the model inclination and error to the plot
    LabelStr="Inc_model \t=\t".expandtabs()+str(AvgModel['INCLINATION'][0])+" $\pm$ " +str(AvgModel['INCLINATION_ERR'][0])+" $^\circ$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    #   Add the global position-angle to the plot along with it's uncertainty.
    LabelStr="PA_model,g \t\t=\t".expandtabs()+str(AvgModel['PA_GLOBAL'])+" $\pm$ " +str(AvgModel['PA_GLOBAL_ERR'])+" $^\circ$"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    #   Add the systemic velocity to the plot
    LabelStr="VSys_model \t=\t".expandtabs()+str(AvgModel['VSYS'][0])+" $\pm$ " +str(AvgModel['VSYS_ERR'][0])+" km/s"
    fig.text(1.2,yTextLoc, LabelStr , ha='left',rotation=0,va='center',size=TextSize,color='k')
    yTextLoc-=YTextStep
    return yTextLoc

    
def AvgPVPlots(fig,ObjDict,AvgModel,CubeInfo,FittingOptions,AnalysisFncs,left,base,w,h,buf):
    """
        This function draws the major and minor axis PV panels onto the diagnostic plot.  It uses routines found in DiagnosticPlots/PVPlotFuncs.py
    """
    #   For the PV plots we need to load in the full resolution MCG model cubelet
    ModelCube=AnalysisFncs['CubeFnc']['BasicCubeAnalysis'](ObjDict['FullResAvgModelCubeName'])

    #   Make the major axis PV plot
    placement=[left+0.1*(w+buf),base+.0*(h+buf),w*0.8,h]
    MajorMinorSwitch=0
    ax=AnalysisFncs['PlotFncs']['AvgModelPVPlot'](fig,placement,ObjDict,AvgModel,CubeInfo,FittingOptions,AnalysisFncs,MajorMinorSwitch,ModelCube)
    #       For the major axis, add in the inclination corrected rotation curve
    AnalysisFncs['PlotFncs']['AddSingleRC_to_PVPlot'](ax,AvgModel,CubeInfo)
    #   Make the minor axis PV plot
    placement=[left+1.1*(w+buf),base+.0*(h+buf),w*0.8,h]
    MajorMinorSwitch=1
    ax=AnalysisFncs['PlotFncs']['AvgModelPVPlot'](fig,placement,ObjDict,AvgModel,CubeInfo,FittingOptions,AnalysisFncs,MajorMinorSwitch,ModelCube)
