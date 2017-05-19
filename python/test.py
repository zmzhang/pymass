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
from sortedcontainers import SortedSet 

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


def find_closest(A, target):
    # A must be sorted
    if len(A) == 0:
        return []
    if len(A) == 1:
        return 0
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A) - 1)
    left = A[idx - 1]
    right = A[idx]
    idx -= target - left < right - target
    return idx.tolist()

def find_idx(rg, rt, mz, threshold):
    rt =rg[find_closest(rg[:,0], rt)][0]
    idx_l = np.searchsorted(rg[:,0], rt, side='left')
    idx_u = np.searchsorted(rg[:,0], rt, side='right')
    idx = find_closest(rg[idx_l:idx_u,1], mz) + idx_l
    if rg[idx,1] > mz-threshold and rg[idx,1] < mz+threshold:
        return idx
    else:
        return -1


def pic_seeds(rmv_sort, mz_tol = 0.5):
    seed_set = SortedSet(key=lambda val: val[1])
    for i in range(rmv_sort.shape[0]):
        idx = seed_set.bisect_left(tuple(rmv_sort[i])) 
        if len(seed_set) == 0:
            seed_set.add(tuple(rmv_sort[0]))
        else:
            if idx == 0:
                if seed_set[0][1] - rmv_sort[i][1] > mz_tol:
                    seed_set.add(tuple(rmv_sort[i]))
            elif idx == len(seed_set):
                if rmv_sort[i][1] - seed_set[-1][1] > mz_tol:
                    seed_set.add(tuple(rmv_sort[i]))
            else:
                if seed_set[idx][1] - rmv_sort[i][1] > mz_tol and \
                   rmv_sort[i][1] - seed_set[idx-1][1] > mz_tol:
                    seed_set.add(tuple(rmv_sort[i]))   
    seeds = np.array(seed_set) 
    return seeds[seeds[:,2].argsort()[::-1],:]

def FPIC(lcms, seed, rt_width, mz_width, b_plot=True):

    pic_ids = []
    stds=[]
    rg   = get_region(seed, lcms, rt_width, mz_width)
    rtm = np.mean(np.diff(lcms.getRT().T))
    pic_ids.append(find_idx(rg, seed[0], seed[1], sys.float_info.max))
    b_left=True
    b_right=True
    for i in range(0,rt_width):
        threshold = sys.float_info.max
        stds.append(np.std(rg[pic_ids,1]))
        if len(stds) >=3:
            #pass
            threshold = 10 * np.mean(stds[1:])
        if b_left:
            rt_left  = rg[pic_ids[0]][0] - rtm
            b_left = find_idx(rg, rt_left, seed[1], threshold) 
            if b_left != -1:
                pic_ids.insert(0, b_left)
            else:
                b_left = False
        
        if b_right:
            rt_right = rg[pic_ids[-1]][0] + rtm 
            idx_r = find_idx(rg, rt_right, seed[1], threshold)
            if idx_r != -1:
                pic_ids.append(idx_r)
            else:
                b_right=False
        
        if b_left == False and b_right == False:
            break;
    pic = np.array([rg[i] for i in pic_ids])
    
    if b_plot:
        plot_region(rg, seed[0], seed[1])  
        plot_region(pic, seed[0], seed[1])
        figure()
        plot(pic[:,0], pic[:,2])
        figure()
        plot(stds)

    return pic

#mzdata2mzxml('F:/resources/MTBLS188/study files/')
#mzfile=u"mixture_bsa300fmol_n3.mzXML"
mzfile=u"MM14_20um.mzxml"

parser=mzXMLParser()


lcms = parser.parseFile(mzfile.encode(sys.getfilesystemencoding()))




rts=lcms.getRT()
bic=lcms.getBIC()
tics=lcms.getTIC()
#plot(rts,bic,'r')
#plot(rts,tics,'g')

rmv = lcms.getAll()
tic()
ids = rmv[:,2].argsort()[::-1]
rmv_sort = rmv[ids,:]
rids = ids.argsort()
toc()

seeds = pic_seeds(rmv_sort[0:10000,:])

pic = FPIC(lcms, seeds[0], 100, 0.5)

print(pic[0][0:3])
print(rmv_sort[rids[int(pic[0][3])]])