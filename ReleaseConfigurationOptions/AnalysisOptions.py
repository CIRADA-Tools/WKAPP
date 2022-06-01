import numpy as np

from . import HydraField as HyaF
from . import HydraDR1_Field as HyaFDR1
from . import NormaField as NorF
from . import NGC4636Field as N4636F

"""
    This module contains routines that control the various fitting options.  These routines point to functions in the individual field modules where the definitions for each field are set.  It contains the routines:
    SetSizeLimit --> Sets the lower limit on Ell_Maj for attempting fitting.
    FATAnalysisOptions --> Sets the options for the FAT runs
    BaroloAnalysisOptions --> Sets the options for the Barolo analysis
    FitComparisonOptions --> Sets the options for the comparing FAT and Barolo runs as well as generating the average models.
"""

def SetSizeLimit():
    """
        Set the lower limit for Ell_Maj in beams for attemping fitting
    """
    SizeLimit=2.
    return SizeLimit

def FATAnalysisOptions(FolderDict):
    """
        Set the options for running FAT --only 1 field should be selected
    """
    #   The Hydra field options
    FATOptions=HyaF.FATAnalysisOptions(FolderDict)
        #   The Norma field options
    #FATOptions=NorF.FATAnalysisOptions(FolderDict)
    #   The NGC4636 field options
    #FATOptions=N4636F.FATAnalysisOptions(FolderDict)
    #   The Hydra DR1 field options
    #FATOptions=HyaFDR1.FATAnalysisOptions(FolderDict)

    return FATOptions

def BaroloAnalysisOptions(FolderDict):
    """
        Set the options for running 3DBarolo --only 1 field should be selected
    """
    #   The Hydra field options
    BaroloOptions=HyaF.BaroloAnalysisOptions(FolderDict)
    #   The Normafield options
    #BaroloOptions=NorF.BaroloAnalysisOptions(FolderDict)
    #   The NGC4636 field options
    #BaroloOptions=N4636F.BaroloAnalysisOptions(FolderDict)
    #   The Hydra DR1 field options
    #BaroloOptions=HyaFDR1.BaroloAnalysisOptions(FolderDict)
    
    return BaroloOptions



def FitComparisonOptions(FolderDict):
    """
       Set the options for comparing FAT and Barolo fits as well as combining the fits into an average model --only 1 field should be selected
    """
 
    #   The Hydra field options
    FittingOptions=HyaF.FitComparisonOptions(FolderDict)
    #   The Norma field options
    #FittingOptions=NorF.FitComparisonOptions(FolderDict)
    #   The NGCC4636 field options
    #FittingOptions=N4636F.FitComparisonOptions(FolderDict)
    #   The Hydra DR1 field options
    #FittingOptions=HyaFDR1.FitComparisonOptions(FolderDict)
    
    return FittingOptions
