# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 22:16:48 2017

@author: zmzhang
"""


import sys
from _pymass import mzXMLParser
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

parser=mzXMLParser()


lcms = parser.parseFile(mzfile.encode(sys.getfilesystemencoding()))




rt=lcms.getRT()
bic=lcms.getBIC()
tics=lcms.getTIC()
plot(rt,bic,'r')
plot(rt,tics,'g')

rg = lcms.getRegion(700,780, 500, 501)

