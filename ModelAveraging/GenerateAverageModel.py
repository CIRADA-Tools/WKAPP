import numpy as np
import astropy
import copy as copy
import os
from scipy import interpolate
import astropy
from astropy.io import fits
from astropy import units as u
from astropy import wcs

from . import CalcSurfDensity as CSD
"""
    This module contains the routines that generate the average model from a set of FAT and 3DBarolo fits.  It contains the routines:
    MakeAvgModel --> This function makes the average model
    SetR --> This function sets the radial grid for the model
    AverageGeometry --> This function averages the geometric parameters
    VelocityAveraging --> This function averages the RC curves
    
    PropagateIncError --> This function propagates the eorrors from the inclination into the rotation curve and surface density profiles
    SetAvgModelConsts --> This function sets a few average model tilted ring parameters to specific constants
    GetAvgModelSkyCoordsAndErrs --> This function converts the pixel centroid coordinates into RA and DEC
    RoundDict --> This function rounds the various tilted ring parameters to a reasonable number of sig figs.
    MakeModelCube --> This function makes an MCG realization of the average model
    GetGlobalPositionAngle --> This function gets the global PA by adjusting the measured PA by the angle the DEC makes with the XY axis.
    AddSoFiAVals --> This function adds the source finding values to the dictionary for cross-matching purposes in later tables
"""

def MakeAvgModel(Fits,FittingOptions,CubeInfo,AnalysisFncs,ObjDict):
    """
        This function averages together a set of different model fits to create an 'average' Tilted Ring model as the WALLABY pilot phase kinematic model
    """
    #   First figure out the radial bins for the averaging
    RUse=SetR(Fits,FittingOptions,CubeInfo)
    #   Next average the geometry
    AvgModel=AverageGeometry(Fits,FittingOptions,RUse)
    #   Now average the inclination corrected and interpolated velocities
    AvgModel=VelocityAveraging(AvgModel,Fits,FittingOptions)
    #   Then get the surface density by ellipse fitting using the average geometry
    AvgModel=CSD.GetNewSurfDens(AvgModel,CubeInfo,AnalysisFncs)
    #   Get the uncertainties due to the inclination
    AvgModel=PropagateIncError(AvgModel)
    #   Set the constant values
    AvgModel=SetAvgModelConsts(AvgModel)
    #   Add some SoFiA values to the average model dictionary for cross-matching in later tables
    AvgModel=AddSoFiAVals(AvgModel,ObjDict)
    return AvgModel
    
def SetR(Fits,FittingOptions,CubeInfo):
    """
        This routine sets the radial grid for the TR parameters.
    """
    #   First make an empty list of minimum and maximum radii
    Rmin=[]
    Rmax=[]
    #   Make a counter for the fit steps
    FitStep=0
    #   Loop through all the fits
    for i in range(FittingOptions['nTotFits']):
        #   Check if the fit was successful
        if Fits[FitStep]['FITAchieved']:
            #   Check if this fit is used for averaging (i.e. only use one of the 2 FAT TR dictionaries made in each FAT analysis).
            if Fits[FitStep]['FitForAveraging']:
                #   Add the smallest and largest radius from each fit to the minimum and maximum lists
                Rmin.append(Fits[FitStep]['R'][0])
                Rmax.append(Fits[FitStep]['R'][-1])
        FitStep+=1
    #   Set the lowest R point as the largest minimum radius of all the fits
    Rlow=np.max(Rmin)
    #   Set the highest R point as the smallest maximum radius of all the fits plus half a beam.  This does a slightly better job of matching and requires only minimal extrapolotions of the rotation curves.
    #       This does require setting the radial bin size as half the beam size
    dR=CubeInfo['CubeHeader']['BMAJ']*3600/2.
    #   Now set Rhigh
    Rhigh=np.min(Rmax)+dR
    #   Make the radial grid
    RArr=np.arange(Rlow,Rhigh,dR)
    #   Return the radial array
    return RArr

def AverageGeometry(Fits,FittingOptions,R):
    """
        This function averages the geometric parameters from a set of different FAT and Barolo fits.  it assumes that the models being averaged are flat disk models
    """
    #   Set the geometry averaging keys
    Keys=['INCLINATION','POSITIONANGLE','XCENTER','YCENTER','VSYS']
    #   Set up an array of summations for each parameter to be averaged.
    Sums=np.zeros(len(Keys))
    #   Set a step counter to figure out where we are in the fitting averages
    FitStep=0
    nSuccess=0
    #   Loop through all the fits
    for i in range(FittingOptions['nTotFits']):
        #   Make sure the fit was successful
        if Fits[FitStep]['FITAchieved']:
            #   If the fit is meant to be used for averaging continue (this avoids double counting FAT fits)
            if Fits[FitStep]['FitForAveraging']:
                #   Loop through all the keywords for averaging
                for j in range(len(Keys)):
                    #   Sum up the geometry value in the 1st radial bin (these should be constant in the fits)
                    Sums[j]+=Fits[FitStep][Keys[j]][0]
                    #   It's possible that the PA may be around 0/360 degrees, making the summation a problem -- check for that
                    if Keys[j] =='POSITIONANGLE':
                        if i >0:
                            #   Get the difference
                            Diff=np.abs(Fits[0][Keys[j]][0]-Fits[FitStep][Keys[j]][0])
                            #   If there is a large difference, we need to fix it
                            if Diff > 180:
                                #   Undo the summation
                                Sums[j]-=Fits[FitStep][Keys[j]][0]
                                #   Adjust the starting PA to be closer to the initial value
                                if Fits[FitStep][Keys[j]][0] >Fits[0][Keys[j]][0]:
                                    Fits[FitStep][Keys[j]][0]-=360
                                else:
                                    Fits[FitStep][Keys[j]][0]+=360.
                                #   Now redo the summation
                                Sums[j]+=Fits[FitStep][Keys[j]][0]
                #   Keep track of the number of fits used in the averaging
                nSuccess+=1
        FitStep+=1
     #       Get the averages
    Avgs=Sums/nSuccess
    #       Now go back through and get the uncertainties
    FitStep=0
    Sums=np.zeros(len(Keys))
        #   Do the same set of checks
    for i in range(FittingOptions['nTotFits']):
        if Fits[FitStep]['FITAchieved']:
            if Fits[FitStep]['FitForAveraging']:
                for j in range(len(Keys)):
                    #   This time look at the difference between the values
                    Sums[j]+=(Fits[FitStep][Keys[j]][0]-Avgs[j])**2.
        FitStep+=1
    #   Get the uncertainties
    Errs=np.sqrt(Sums/(nSuccess))
    #   Due to the possible PA issues, make sure the sum is between 0-360
    if Avgs[1] < 0:
        Avgs[1]+=360.
    elif Avgs[1] > 360.:
        Avgs[1]-=360.
    #   Put all the averages into a Tilted Ring dictionary
    #       Start with an empty model dictionary
    AvgGeometryDict={}
    #       Loop through all the key words
    for i in range(len(Keys)):
        #   Extend the average across all model radial bins using the R/R notation
        AvgGeometryDict[Keys[i]]=R/R*Avgs[i]
        #   Do the same thing for the uncertainties
        ErrKey=Keys[i]+'_ERR'
        AvgGeometryDict[ErrKey]=R/R*Errs[i]
    #   Set the average model radius and surface density radial grids
    AvgGeometryDict['R']=R
    AvgGeometryDict['R_SD']=R
    #   Return the average model dictionary
    return AvgGeometryDict

def VelocityAveraging(AvgModel,Fits,FittingOptions):
    """
        This function averages the inclination corrected and interpolated velocity profiles into a model rotation curve
    """
    #   Set up the starting parts of the loop
    FitStep=0
    nSuccess=0
    #   Initialize an empty array for intepolated velocities
    VInterpolate=[]
    #   Loop through all models
    for i in range(FittingOptions['nTotFits']):
        #   Check that the fit was successful
        if Fits[FitStep]['FITAchieved']:
            #   Check if the fit should be used for averaging (this stops double counting FAT fits)
            if Fits[FitStep]['FitForAveraging']:
                #   Place VRot into a working array
                VRot_Work=copy.copy(Fits[FitStep]['VROT'])
                #   Adjust VRot to the new average inclination
                VRot_Work=VRot_Work/np.sin(AvgModel['INCLINATION'][0]*np.pi/180.)*np.sin(Fits[FitStep]['INCLINATION'][0]*np.pi/180.)
                #   Use the scipy interpolate term to generate an interpolation function
                V_IntFunc=interpolate.UnivariateSpline(Fits[FitStep]['R'], VRot_Work,k=2)
                #   Store the interpolated array for calculation of the uncertainty
                #       This call makes VInterpolate a 2D array with the shape (nSuccess, nR)
                VInterpolate.append(V_IntFunc(AvgModel['R']))
                #   Add the inclination corrected intepolated rotation curve to a VAvg array
                if nSuccess == 0:
                    VAvg=copy.copy(VInterpolate[nSuccess])
                else:
                    VAvg+=VInterpolate[nSuccess]
                #   Keep track of the number of fits that go into the averaging
                nSuccess+=1
        FitStep+=1
    #   Average the profiles
    VAvg=VAvg/nSuccess
    #   Now get the uncertainties by looping through all the successful profiles stored in the VInterpolate array
    for i in range(nSuccess):
        if i == 0:
            VErr=(VInterpolate[i]-VAvg)**2.
        else:
            VErr+=(VInterpolate[i]-VAvg)**2.
    #   Set the uncertainty
    VErr=np.sqrt(VErr/nSuccess)
    #   Save the average velocities and errors to the tilted ring model
    AvgModel['VROT']=VAvg
    AvgModel['VROT_ERR']=VErr
    #   Return the tilted ring model
    return AvgModel

def PropagateIncError(AvgDict):
    """
        This function progates errors in the inclination to errors in the rotation and surface density profiles
    """
    #   Start with the model inclination
    Inc=AvgDict['INCLINATION']*np.pi/180.
    #   Get the smallest and largest inclinations
    IncMin=(AvgDict['INCLINATION']-AvgDict['INCLINATION_ERR'])*np.pi/180.
    IncMax=(AvgDict['INCLINATION']+AvgDict['INCLINATION_ERR'])*np.pi/180.
    #   Get the rotation curves for the models with the smaller and larger inclinations
    VRotMin=copy.copy(AvgDict['VROT'])*np.sin(Inc)/np.sin(IncMin)
    VRotMax=copy.copy(AvgDict['VROT'])*np.sin(Inc)/np.sin(IncMax)
    #   Get the difference between the minimum and maximum rotations
    DeltaV=abs(VRotMax-VRotMin)
    #   Set the error to half the difference between the two measurements
    AvgDict['VROT_INC_ERR']=DeltaV/2.
    #   Do the same thing for the surface density
    SDMin=copy.copy(AvgDict['SURFDENS_FACEON'])*np.cos(IncMin[0])/np.cos(Inc[0])
    SDMax=copy.copy(AvgDict['SURFDENS_FACEON'])*np.cos(IncMax[0])/np.cos(Inc[0])
    DeltaSD=abs(SDMax-SDMin)
    AvgDict['SURFDENS_INC_ERR']=DeltaSD/2.
    return AvgDict

def SetAvgModelConsts(AvgModel):
    """
        This function sets some of the tilted ring dictionary parameters that are not fit by either the FAT or Barolo runs.  The radial velocity is set to 0, the height to 0, and the dispersion to 10 km/s
    """
    AvgModel['VRAD']=0.*AvgModel['R']/AvgModel['R']
    AvgModel['Z0']=0.*AvgModel['R']/AvgModel['R']
    AvgModel['VDISPERSION']=10.*AvgModel['R']/AvgModel['R']
    return AvgModel


def GetAvgModelSkyCoordsAndErrs(AvgModel,CubeInfo,AstroFncs):
    """
        This function gets the average model's sky coordinates and uncertainties using the centroid pixels, the cube header (in principle, this should be the original SoFiA cubelet to avoid the fisheye effect) and astropy wcslib routines
    """
    #   This requires the center point to be set
    X=AvgModel['XCENTER']
    Y=AvgModel['YCENTER']
    #   Get the low and high values for X and Y based on the uncertainties.
    XLow=X[0]-AvgModel['XCENTER_ERR'][0]
    XHigh=X[0]+AvgModel['XCENTER_ERR'][0]
    YLow=Y[0]-AvgModel['YCENTER_ERR'][0]
    YHigh=Y[0]+AvgModel['YCENTER_ERR'][0]
    #   Set up definitions for the low and high coordinates
    CoordLow=[[XLow,YLow,0]]
    CoordHigh=[[XHigh,YHigh,0]]
    #   Use wcslib to convert these pixels into RA and DEC
    PosLow=CubeInfo['CubeWCS'].wcs_pix2world(CoordLow,0)
    PosHigh=CubeInfo['CubeWCS'].wcs_pix2world(CoordHigh,0)
    #   Get RA and DEC for the X and Y positions (see BasicStatistics/AstroFncs.py)
    AvgModel['RA'],AvgModel['DEC']=AstroFncs['CalcRA_Dec_FromCube'](X,Y,CubeInfo)
    #   Set the uncertainties as half the difference between the low and high measures
    AvgModel['RA_ERR']=abs(PosLow[0,0]-PosHigh[0,0])/2.
    AvgModel['DEC_ERR']=abs(PosLow[0,1]-PosHigh[0,1])/2.
    return AvgModel

def RoundDict(AvgDict):
    """
        This routine is meant to round all the dictionary outputs to reasonable numbers of sig-figs
    """
#   Loop through all radii
    for i in range(len(AvgDict['R'])):
        #   Round all radial profiles appropriately
        AvgDict['R'][i]=round(AvgDict['R'][i],2)
        AvgDict['INCLINATION'][i]=round(AvgDict['INCLINATION'][i],1)
        AvgDict['INCLINATION_ERR'][i]=round(AvgDict['INCLINATION_ERR'][i],1)
        AvgDict['POSITIONANGLE'][i]=round(AvgDict['POSITIONANGLE'][i],1)
        AvgDict['POSITIONANGLE_ERR'][i]=round(AvgDict['POSITIONANGLE_ERR'][i],1)
        AvgDict['VSYS'][i]=round(AvgDict['VSYS'][i],1)
        AvgDict['VSYS_ERR'][i]=round(AvgDict['VSYS_ERR'][i],1)
        AvgDict['XCENTER'][i]=round(AvgDict['XCENTER'][i],3)
        AvgDict['XCENTER_ERR'][i]=round(AvgDict['XCENTER_ERR'][i],4)
        AvgDict['YCENTER'][i]=round(AvgDict['YCENTER'][i],3)
        AvgDict['YCENTER_ERR'][i]=round(AvgDict['YCENTER_ERR'][i],4)
        AvgDict['VROT'][i]=round(AvgDict['VROT'][i],1)
        AvgDict['VROT_ERR'][i]=round(AvgDict['VROT_ERR'][i],1)
        AvgDict['VROT_INC_ERR'][i]=round(AvgDict['VROT_INC_ERR'][i],1)
        #       For the RA and Dec, round to 4 decimal placess
        AvgDict['RA'][i]=round(AvgDict['RA'][i],6)
        AvgDict['DEC'][i]=round(AvgDict['DEC'][i],6)
    #   Do the extra rounding needed for the surface density profiles as they are on a separate radial grid
    for i in range(len(AvgDict['R_SD'])):
        #   For the surface density use 2 decimal places
        AvgDict['R_SD'][i]=round(AvgDict['R_SD'][i],2)
        AvgDict['SURFDENS'][i]=round(AvgDict['SURFDENS'][i],2)
        AvgDict['SURFDENS_ERR'][i]=round(AvgDict['SURFDENS_ERR'][i],2)
        AvgDict['SURFDENS_INC_ERR'][i]=round(AvgDict['SURFDENS_INC_ERR'][i],3)
        AvgDict['SURFDENS_FACEON'][i]=round(AvgDict['SURFDENS_FACEON'][i],3)
    #       Note that RA_Err and Dec_Err are scalars so place them out of the loop
    AvgDict['RA_ERR']=round(AvgDict['RA_ERR'],6)
    AvgDict['DEC_ERR']=round(AvgDict['DEC_ERR'],6)
    AvgDict['PA_GLOBAL']=round(AvgDict['PA_GLOBAL'],1)
    AvgDict['PA_GLOBAL_ERR']=round(AvgDict['PA_GLOBAL_ERR'],1)
    return AvgDict
    
def MakeModelCube(AvgModel,MCGFncs,CubeFncs,ObjDict,SmoothSwitch):
    """
        This function makes an MCG cubelet realization of the tilted ring model -- see MCGModel/MCGModelGeneration.py for some of the calls
    """
    #   Set the MCG dictionary names
    MCGDict=MCGFncs['MCG_Defs'](ObjDict['ObjName_Underscored'])
    #   Select if we are making a lower res or higher res cube
    if SmoothSwitch == 0:
        TargCubeName=ObjDict['VelocityCubeFileName']
    elif SmoothSwitch == 1:
        TargCubeName=ObjDict['SmoothedVelocityCubeFileName']
    #   import astropy fits routines
    import astropy
    from astropy.io import fits
    from astropy import units as u
    from astropy import wcs
    #   Load in the needed cube header
    Cube=fits.open(TargCubeName)
    CubeHeader=Cube[0].header
    CubeWCS = wcs.WCS(CubeHeader)
    #   Get the start location of the pixels from the header
    CubeHeader=CubeFncs['GetStartPixelLoc'](CubeWCS,CubeHeader)
    #   Close the cube
    Cube.close()
    #   Make the MCG Cube
    MCGFncs['MakeMCGModel'](MCGDict,AvgModel,CubeHeader)
    return MCGDict

def CompareAvgModelToFitData(AvgModel,ObjDict,AvgOutputDict,FittingOptions,CubeFnc):
    """
        This function compares a model cube to the underlying data by generating a comparison cube.
    """
    #   Name the difference cube
        #   Set the base name for output files
    KinTR_Underscore=FittingOptions['ProvenanceVals']['KINTR'].split()
    KinTR_Underscore=KinTR_Underscore[0]+"_"+KinTR_Underscore[1]+"_"+KinTR_Underscore[2]
    OutputBaseName=ObjDict['ObjName_Underscored']+"_"+KinTR_Underscore
    #       Nose set the output name
    DiffCubeName=AvgOutputDict['OutputDirName']+OutputBaseName+"_DifferenceCube.fits"
    #   Load in the smoothed cubelet and the data cube
    SmoothDataCube=CubeFnc['BasicCubeAnalysis'](ObjDict['SmoothedVelocityCubeFileName'])
    #   Open the model cube with the same resolution
    ModelCube=fits.open(ObjDict['SmoothedResAvgModelCubeName'])
    #   Calculate the difference between the model and the data
    ModelCube[0].data=ModelCube[0].data-SmoothDataCube['Data']
    #   Write out the difference cube
    ModelCube.writeto(DiffCubeName,overwrite=True)
    #   Close the model cube
    ModelCube.close()
    
def GetGlobalPositionAngle(AvgModel,CubeInfo,AstroFncs):
    """
       The kinematic modelling codes get the position angle relative to the pixel axes, but these may be tilted with respect to RA and DEC.  This routine gets the position angle in the RA and DEC coordinates using the astropy wcslib functions.
    """
    #   Set the location of the center in pixels as well as the average model RA and DEC
    XPix=AvgModel['XCENTER'][0]
    YPix=AvgModel['YCENTER'][0]
    RA_Cent=AvgModel['RA'][0]
    DEC_Cent=AvgModel['DEC'][0]
    #   Set a new DEC value one beam higher
    DEC_Delt=DEC_Cent+CubeInfo['CubeHeader']['BMAJ']
    #   Set the RA_Cent and DEC_Delt into a coordinate
    CoordDelt=[[RA_Cent,DEC_Delt,0]]
    #   Convert this point to pixel coordinates.
    PosDelt=CubeInfo['CubeWCS'].wcs_world2pix(CoordDelt,0)
    #   Now figure out how far away this new point is from the center point
    XDelt=PosDelt[0][0]-XPix
    YDelt=PosDelt[0][1]-YPix
    #   Use these points to get the angle the declination makes with the X-Y axis
    DEC_Angle=np.arctan2(YDelt,XDelt)*180./np.pi
    #   Get the global PA by adjusting the XY version of PA by the angle DEC makes with the XY axis
    PA_Global=AvgModel['POSITIONANGLE'][0]-(DEC_Angle-90.)
    #   Ensure that the PA is between 0 and 360
    if PA_Global < 0.:
        PA_Global+=360.
    elif PA_Global > 360.:
        PA_Global-=360.
    #   Store the global PA and associated error into the tilted ring dictionary
    AvgModel['PA_GLOBAL']=PA_Global
    AvgModel['PA_GLOBAL_ERR']=AvgModel['POSITIONANGLE_ERR'][0]
    return AvgModel

def AddSoFiAVals(AvgModel,ObjDict):
    """
        This function copies some of the SoFiA values into the average model dictionary so that they can be stored in a later table
    """
    AvgModel['ra_sofia']=ObjDict['CatEntry']['ra']
    AvgModel['dec_sofia']=ObjDict['CatEntry']['dec']
    AvgModel['freq_sofia']=ObjDict['CatEntry']['freq']
    return AvgModel
