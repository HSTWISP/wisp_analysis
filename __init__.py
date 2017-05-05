import copy as cp
import scipy.signal as si 
import numpy as np
import matplotlib.pyplot as plt
import pickle

import os
import astropy.io.ascii as asciitable
import astropy.io.fits as fits 
#import time  
import mpfit
#import multiprocessing as mp
import math
from scipy import interpolate
from scipy import integrate

from trim_spec import trim_spec
from utilities import gaussian 
from utilities import is_number
from utilities import read_config
#from specmodel import emissionline_model_spline 
#from specmodel import model_resid_spline 
from find_cwt import find_cwt  
from find_cwt  import loop_field_cwt
from find_cwt import test_obj_cwt 
from fitting import emissionline_model
from fitting import model_resid 
from fitting import fit_obj
from fitting import fitandplot
from guis import * 
from measure_z_interactive import * 
import pickle

try:
    from stacking import * 
except ImportError:
    pass
#    print 'No stacking module. It is not needed for line finding. Skipping'

