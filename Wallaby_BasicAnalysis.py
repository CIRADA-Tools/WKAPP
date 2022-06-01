import numpy as np
import pandas as pd

import ReleaseConfigurationOptions as RCO
import SimpleIO as SIO
import BasicStatistics as BS
import CubeLoad as CL

"""
    This script does a basic set of measurements on each galaxy in a field (size in beams, mass, predicted R_HI, # of velocity elements, volume estimate, and S/N).  Some of these are used for the Barolo and FAT analysis selections in later steps in the proto-pipeline.
"""

def Main():
    #   Get the various file and folder names
    FolderDict=RCO.FF.BasicFolderAndFileLoc()
    #   Get the details about the excel file
    ExcelDetailsDict=RCO.FF.ExcelCatDetails()
    #   Load in the general WALLABY catalogue
    WallabyCat=SIO.SCIO.LoadXLS(FolderDict['CatalogueFile'],ExcelDetailsDict)
    #   Get the number of galaxies
    nGalaxies=np.shape(WallabyCat['id'])[0]
    #   Convert the ID column into integers
    WallabyCat['id']=WallabyCat['id'].astype('int')
    #   Add the frequency bin size to the catalogue
    WallabyCat['df']=1.85185185185E+04
 
    #   Get the object sizes in beams using the SoFiA Ell_Maj parameter
    WallabyCat['ObjectSize_Beams']=BS.BOM.SizeEstimate(WallabyCat)
    #   Get the distance to the object (assuming H0=70)
    WallabyCat['Distance']=BS.BOM.DistanceEstimate(WallabyCat)
    #   Get the HI mass in the object
    WallabyCat['logM']=BS.BOM.MassCalc(WallabyCat)
    #   Get the predicted size of the object
    WallabyCat['RHI'],WallabyCat['RHI_Projected']=BS.BOM.RHI_Estimate(WallabyCat)
    #   Get the number of velocity channels
    WallabyCat['VelElements']=BS.BOM.VelElements(WallabyCat)
    #   Get an estimate of the volume of the object
    WallabyCat['VolEstimate']=BS.BOM.CubeVolEstimate(WallabyCat)

    #   Set the up the nosie arrays
    RMS=np.zeros(nGalaxies)
    SN_Peak=np.zeros(nGalaxies)
    SN_Avg=np.zeros(nGalaxies)
    SN_Obs=np.zeros(nGalaxies)

    #   Get the S/N ratio for each object
    
    #for i in range(1,2):
    for i in range(1,nGalaxies+1):
        #   Get the basic object information
        j=i-1
        #   Set up the specific galaxy dictionary object by the step
        ObjDict=RCO.FF.InitializeObjectDictionary(WallabyCat,i,FolderDict)
        #   Estimate the S/N of the cube in a number of different ways
        RMS[j],SN_Peak[j],SN_Avg[j],SN_Obs[j]= BS.BOM.SNEstimate(ObjDict,WallabyCat)

    #   Place everything we want to save into a large dictionary
    MeasurementDict={'NAME':WallabyCat['name'],'ID':WallabyCat['id'],'SIZE':WallabyCat['ObjectSize_Beams'],'DISTANCE':WallabyCat['Distance'],'HI_FLUX':WallabyCat['f_sum'],'LOG_MASS':WallabyCat['logM'],'R_HI':WallabyCat['RHI'],'R_HI_PROJECTED':WallabyCat['RHI_Projected'],'VEL_ELEMENTS':WallabyCat['VelElements'],'VOLUME_ESTIMATE':WallabyCat['VolEstimate'],'RMS':RMS,'SN_PEAK':SN_Peak,'SN_AVG':SN_Avg,'SN_OBS':SN_Obs}
  
    #   Save the dictionary to file
    DF=pd.DataFrame.from_dict(MeasurementDict)
    csvname=FolderDict['MeasurementCatalogue']
    DF.to_csv(csvname,index=False)

Main()
