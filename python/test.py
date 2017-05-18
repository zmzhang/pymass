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


def plot_region(rmv, lcms, n):
    (rt, mz, val) = rmv[n]
    region = lcms.getRegion(rt - 200, rt + 200, mz - 0.5, mz + 0.5)
    rg = np.array(region).reshape((len(region), region[0].shape[0]))
    figure()
    scatter(rt, mz, c = 'r', marker = 'x', s = 100)
    scatter(rg[:,0], rg[:,1], c = np.log(rg[:,2]))
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

plot_region(rmv_sort, lcms, 0)



(rt, mz, val) = rmv[2]
region = lcms.getRegion(rt - 100, rt + 100, mz - 0.3, mz + 0.3)

