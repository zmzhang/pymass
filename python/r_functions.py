# -*- coding: utf-8 -*-
"""
Created on Mon May 29 10:46:39 2017

@author: humbl
"""
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
import os
import pandas as pd
import numpy as np

def init_r():
    """
        1. Install R 3.4.1 with conda
        2. Install XCMS and PITracer in R
        3. Install rpy2 in with conda
    """
    numpy2ri.activate()
    wd = os.path.dirname(os.path.abspath(__file__)).replace('\\','/')
    robjects.r(f'''setwd(\'{wd}\')''')
    robjects.r('''source('r_functions.R')''')
    PIT = robjects.globalenv['PIT']
    XC = robjects.globalenv['XC']
    return PIT, XC

def XCMS(mzMLFile, w1, w2, snr, intensity):
    df = pd.DataFrame(columns=['rt', 'mz', 'intensity'])
    _,XC = init_r()
    peaks_xcms = np.array(XC(mzXMLFile=mzMLFile, w1=w1, w2=w2, snr=snr, intensity=intensity))
    for i in range(peaks_xcms.shape[0]):
        rt  = peaks_xcms[i,3]
        mz  = peaks_xcms[i,0]
        intensity = peaks_xcms[i,6]
        df.loc[len(df)] = [rt, mz, intensity] 
    return df