import numpy as np
import pandas as pd
import os.path
from os import path

"""
    This module contains routines for reading in CSV and Excel catalogue files.  It has the routines:
        LoadCSV --> This function loads in a CSV file and puts the catalogue into a dictionary
        LoadXLS --> This function loads in an excel file and puts it into a dictionary
"""

def LoadCSV(File,sep):
    """
        This function loads in a CSV file and puts the catalogue into a dictionary
    """
    #   Read the CSV file
    Contents=pd.read_csv(File,sep=sep)
    #   Set the file contents to a dictionary
    Contents = Contents.rename(columns=lambda x: x.strip())
    return Contents
    
def LoadXLS(File,DetailsDict):
    """
        This function loads in an excel file and puts it into a dictionary
    """
    #   Read in the excel file
    xls = pd.read_excel(File,sheet_name=DetailsDict['SheetName'],header=DetailsDict['HeaderLen'])
    #   If there are rows that need to be deleted due to nan's, give the
    #       column name to use for that stripping
    if DetailsDict['StripCol'][0]== 1:
        xls=xls[xls[DetailsDict['StripCol'][1]].notna()]
    #   Set the file contents to a dictionary
    xls = xls.rename(columns=lambda x: x.strip())
    return xls
