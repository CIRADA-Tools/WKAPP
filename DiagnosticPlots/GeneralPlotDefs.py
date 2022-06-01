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

"""
    This module contains a number of plotting functions that are useful for both diagnostic and average model plots.  This module contains more general definitions.  It contains the routines:
    ColorSelection --> This function automatically selects the line colors for a set of models.
"""
        
def ColorSelection(LineNum):
    """
        This function selects the line color for some model -- it does this by a number.
    """
    if LineNum == -1:
        col='black'
    elif LineNum == 0:
        col='red'
    elif LineNum == 1:
        col='blue'
    elif LineNum == 2:
        col='green'
    elif LineNum == 3:
        col='magenta'
    elif LineNum == 4:
        col='purple'
    elif LineNum == 5:
        col='orange'
    elif LineNum == 6:
        col='brown'
    elif LineNum == 7:
        col='deeppink'
    elif LineNum == 8:
        col='olive'
    elif LineNum == 9:
        col='royalblue'
    elif LineNum == 10:
        col='indigo'
    elif LineNum == 11:
        col='salmon'
    elif LineNum == 12:
        col='teal'
    elif LineNum == 13:
        col='forestgreen'
    return col

