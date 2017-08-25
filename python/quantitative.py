# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 09:39:56 2017

@author: zmzhang
"""

import pandas as pd
import subprocess, os, sys
import pyopenms
from _pymass import mzXMLParser
import _pymass as pm
from FPIC import data2mzxml, tic, toc
import numpy as np

df_std = pd.DataFrame(columns=['name','rt', 'mz']) 
df_std.loc[len(df_std)] = ['propionyl', 1.4, 221.15751] 
df_std.loc[len(df_std)] = ['nialamide', 5.7, 299.15025] 
df_std.loc[len(df_std)] = ['sulfadimethoxine', 8.4, 317.11851]
df_std.loc[len(df_std)] = ['reserpine', 10.6, 609.28065]
df_std.loc[len(df_std)] = ['terfenadine', 12.4, 472.32099]
df_std.loc[len(df_std)] = ['hexadecanoyl', 16.3, 403.36096]
df_std.loc[len(df_std)] = ['octadecanoyl', 18.1, 431.39227]

df_std.rt = df_std.rt * 60


def FeatureFindingMetabo(mzfile, noise_threshold_int, snr):
    finder = 'C:/Program Files/OpenMS/bin/FeatureFinderMetabo.exe'
    feature_file = 'tmp.featureXML'
    noise_threshold_int = noise_threshold_int / snr
    subprocess.call([finder, '-in', mzfile, '-out', feature_file, 
               '-algorithm:common:noise_threshold_int', f'{noise_threshold_int}',
               '-algorithm:common:chrom_peak_snr', f'{snr}',
               '-algorithm:common:chrom_fwhm', '10',
               '-algorithm:mtd:mass_error_ppm', '20',
               '-algorithm:mtd:reestimate_mt_sd', 'true',
               '-algorithm:mtd:min_sample_rate', '0',
               '-algorithm:mtd:min_trace_length', '2',
               '-algorithm:epd:width_filtering', 'off',
               '-algorithm:ffm:charge_lower_bound', '1',
               '-algorithm:ffm:charge_lower_bound', '5'])  
    featuremap = pyopenms.FeatureMap()
    featurexml = pyopenms.FeatureXMLFile()
    featurexml.load(feature_file, featuremap)
    os.remove(feature_file)
    return featuremap


def parse_featureXML_FFM(featuremap):   
    df = pd.DataFrame(columns=['rt', 'mz', 'intensity'])   
    for i in range(featuremap.size()):
        feature = featuremap[i]
        isotope_distances = feature.getMetaValue(b'isotope_distances')
        rt = feature.getRT()
        mz = feature.getMZ()
        intensity = feature.getIntensity()
        for j in range(feature.getMetaValue(b'num_of_masstraces')):
            if j == 0:
                df.loc[len(df)] = [rt, mz, intensity]
            else:
                mz_delta = isotope_distances[j-1]
                mz = mz + mz_delta
                df.loc[len(df)] = [rt, mz, intensity] 
    return df

def pics2df(pics):
    df = pd.DataFrame(columns=['rt', 'mz', 'intensity'])
    for i,pic in enumerate(pics):
        idx = pic[:,2].argmax()
        rt  = pic[idx,0]
        mz  = pic[idx,1]
        intensity = pic[idx,2]
        df.loc[len(df)] = [rt, mz, intensity] 
    return df

cs = ['1000', '0500', '0200', '0100', '0050', '0020', '0010', '0005', '0002', '0001', '0000']

result_ffm = pd.DataFrame(np.float64('NaN'),index = range(len(cs)),
                          columns=df_std.name.tolist() + ['concentration']) 
result_fpic = pd.DataFrame(np.float64('NaN'),index = range(len(cs)), 
                           columns=df_std.name.tolist() + ['concentration']) 

for k,c in enumerate(cs):
    result_ffm.at[k,'concentration'] = c
    result_fpic.at[k,'concentration'] = c
    name = f'2012_02_03_PStd_{c}_3'
    mzMLfile = f'MTBLS234/{name}.mzML'
    mzXMLfile = f'MTBLS234/{name}.mzxml'
    intensity = 30   
    feature_map = FeatureFindingMetabo(mzMLfile, intensity, 3)
    df_ffm = parse_featureXML_FFM(feature_map)
    
    tic()
    parser=mzXMLParser()
    lcms = parser.parseFile(mzXMLfile.encode(sys.getfilesystemencoding()))
    pics_c = pm.FPICs(lcms, intensity, 200.0, 0.5)
    toc()
    df_fpic = pics2df(pics_c)
    
    for i in range(len(df_std)):
        mz = df_std.at[i, 'mz']
        rt = df_std.at[i, 'rt']
        mz_tol = 0.014
        if i == 6: # peak shape of Octadecanoyl-L-carnitine-d3 is bad
            rt_tol = 60
        else:
            rt_tol = 30
        mz_l = mz - mz_tol
        mz_u = mz + mz_tol
        rt_l = rt - rt_tol
        rt_u = rt + rt_tol
        query_str = f'{mz_l} < mz < {mz_u} and {rt_l} < rt < {rt_u}'
        
                      
        df_ffm_filtered = df_ffm.query(query_str)
        if len(df_ffm_filtered) >0:
            result_ffm.at[k, df_std.at[i,'name']] = df_ffm_filtered.intensity.values.max()

        df_fpic_filtered = df_fpic.query(query_str)
        if len(df_fpic_filtered) >0:
            result_fpic.at[k, df_std.at[i,'name']] = df_fpic_filtered.intensity.values.max()

    

cc_ffm_rd = result_ffm.corr().round(4)
cc_fpic_rd = result_fpic.corr().round(4)
result_ffm_rd =  result_ffm.round(0)
result_fpic_rd = result_fpic.round(0)
