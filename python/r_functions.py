# -*- coding: utf-8 -*-
"""
Created on Mon May 29 10:46:39 2017

@author: humbl
"""
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
import os

def init_r():
    numpy2ri.activate()
    wd = os.path.dirname(os.path.abspath(__file__)).replace('\\','/')
    robjects.r(f'''setwd(\'{wd}\')''')
    robjects.r('''source('r_functions.R')''')
    PITracer = robjects.globalenv['PIT']
    XCMS = robjects.globalenv['XC']
    return PITracer, XCMS

