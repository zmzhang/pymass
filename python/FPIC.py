# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 22:16:48 2017

@author: zmzhang
"""


import sys
from _pymass import mzXMLParser
import _pymass as pm
import numpy as np
from pylab import plot, show, figure, scatter, xlabel, ylabel, hist
import pylab
from matplotlib.ticker import FormatStrFormatter
import subprocess, os

def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")


def get_region(seed, lcms, width, mz_tol):
    (rt, mz, val) = seed
    region = lcms.getRegion(rt - width, rt + width, mz - mz_tol, mz + mz_tol)
    return np.array(region).reshape((len(region), region[0].shape[0])) 

def plot_region(rg, rt = None, mz = None):
    figure()
    scatter(rg[:,0], rg[:,1], c = np.log(rg[:,2]), s = 1)
    if rt != None and mz != None:
        scatter(rt, mz, c = 'r', marker = 'x', s = 100)
    xlabel('Retention Time(S)')
    ylabel('M/Z')
    ax = pylab.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    show()

def mzdata2mzxml(path, converter = 'C:/Program Files/OpenMS/bin/FileConverter.exe'):
    files=os.listdir(path)
    for f in files:
        if f.lower().endswith(".mzdata"): 
            file_in  = path + f
            file_out = path + f[0:-6] + "mzxml"
            subprocess.Popen([converter, '-in', file_in, '-out', file_out])


def pics_id(pics, n):
    pic_ids=[]
    for pic in pics:
        pic_ids.append(pic[pic[:,2].argmax(), 3])
    pic_ids_a = np.array(pic_ids)
    pic_ids_sort = pic_ids_a[pic_ids_a.argsort()].astype(np.int).tolist()
    if len(pics)<n:
        pic_ids_sort = pic_ids_sort + [0]*(n-len(pics))
    return pic_ids_sort

def pics2peaks(pics):
    peaks = np.zeros((len(pics), 7))
    for i, pic in enumerate(pics):
        idx = pic[:,2].argmax()
        peaks[i,0] = pic[idx, 1]         # MZ apex
        peaks[i,1] = np.min(pic[:, 1])   # MZ min
        peaks[i,2] = np.max(pic[:, 1])   # MZ max
        peaks[i,3] = pic[idx, 0]         # RT apex
        peaks[i,4] = pic[0  , 0]         # RT min
        peaks[i,5] = pic[-1 , 0]         # RT max
        peaks[i,6] = pic[idx, 2]         # Intensity max   
    return peaks

if __name__=="__main__":
    
    #mzdata2mzxml('F:/resources/MTBLS188/study files/')
    mzfile=u"MM14_20um.mzxml"
    parser=mzXMLParser()
    lcms = parser.parseFile(mzfile.encode(sys.getfilesystemencoding()))
    pics_c = pm.FPICs(lcms, 300.0, 100.0, 0.5)
    
    
    
    figure()
    rmv      = lcms.getAll()
    ids      = rmv[:,2].argsort()[::-1]
    rmv_sort = rmv[ids,:]
    hist(np.log10(rmv_sort[:,2]), bins = 800)
    
    
    seed = rmv_sort[0,:]
    stds = pm.FPICStd(lcms, seed, 100, 0.5)
    pic = pm.FPIC(lcms, seed, 100, 0.5)
    rg =  get_region(seed, lcms, 100, 0.5)
    plot_region(rg)
    figure()
    plot(pic[:,0], pic[:,2])
    figure()
    plot(stds)