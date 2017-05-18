# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 22:16:48 2017

@author: zmzhang
"""


import sys
from _pymass import mzXMLParser
import numpy as np
from pylab import plot, show, figure, scatter, xlabel, ylabel
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


def get_region(rmv, lcms, n, width, mz_tol):
    (rt, mz, val) = rmv[n]
    region = lcms.getRegion(rt - width, rt + width, mz - mz_tol, mz + mz_tol)
    return np.array(region).reshape((len(region), region[0].shape[0])) 

def plot_region(rg, rt = None, mz = None):
    figure()
    if rt != None and mz != None:
        scatter(rt, mz, c = 'r', marker = 'x', s = 100)
    scatter(rg[:,0], rg[:,1], c = np.log(rg[:,2]), s = 1)
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


#mzdata2mzxml('F:/resources/MTBLS188/study files/')


mzfile=u"mixture_bsa300fmol_n3.mzXML"

parser=mzXMLParser()


lcms = parser.parseFile(mzfile.encode(sys.getfilesystemencoding()))




rt=lcms.getRT()
bic=lcms.getBIC()
tics=lcms.getTIC()
plot(rt,bic,'r')
plot(rt,tics,'g')

rmv = lcms.getAll()
tic()
rmv_sort = rmv[rmv[:,2].argsort()[::-1],:]
toc()

i = 0
rg = get_region(rmv_sort, lcms, i, 100, 0.5)
plot_region(rg, rmv_sort[i][0], rmv_sort[i][1])

