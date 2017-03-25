# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 22:16:48 2017

@author: zmzhang
"""


import sys
from _pymass import MZXML
from matplotlib.pylab import plot, show, figure


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


mzfile=u"mixture_bsa300fmol_n3.mzXML"

mz1=MZXML()

tic()
mz1.parseFile(mzfile.encode(sys.getfilesystemencoding()))
toc()



rt=mz1.getRT()
bic=mz1.getBIC()
tic=mz1.getTIC()
plot(rt,bic,'r')
plot(rt,tic,'g')

mz=mz1.getMS(0)
val=mz1.getVal(0)

figure()
plot(mz,val)
show()
