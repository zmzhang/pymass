# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 08:14:15 2017

@author: zmzhang
"""


import subprocess, sys, os
import pyopenms
import pandas as pd
from FPIC import data2mzxml, tic, toc
from _pymass import mzXMLParser
import _pymass as pm
import numpy as np


def simulation(fasta, contaminants, out, out_cntm, stddev,
               simulator = 'C:/Program Files/OpenMS/bin/MSSimulator.exe'):   
    """
        Should copy "C:\Program Files\OpenMS\share\OpenMS\examples" to working directory of Python
    """ 
   
    subprocess.call([simulator, '-in', fasta, '-out', out, '-out_cntm',out_cntm, 
               '-algorithm:MSSim:RawSignal:contaminants:file', contaminants,
               '-algorithm:MSSim:RawSignal:noise:detector:stddev', f'{stddev}',
               '-algorithm:MSSim:RawSignal:resolution:value', '5000',
               '-algorithm:MSSim:RawSignal:resolution:type', 'constant',
               '-algorithm:MSSim:Ionization:mz:lower_measurement_limit', '10',
               '-algorithm:MSSim:Ionization:mz:upper_measurement_limit', '1000',
               '-algorithm:MSSim:RT:total_gradient_time', '1000',
               '-algorithm:MSSim:RT:sampling_rate', '0.25',
               '-algorithm:MSSim:RT:scan_window:min', '0',
               '-algorithm:MSSim:RT:scan_window:max', '1000'])

def parse_featureXML_GT(feature_file):
    featuremap = pyopenms.FeatureMap()
    featurexml = pyopenms.FeatureXMLFile()
    featurexml.load(feature_file, featuremap)
    
    hulls = pd.DataFrame(columns=['rt_min', 'rt_max', 'mz_min', 'mz_max', 'detected', 'pic_id'])   
    for i in range(featuremap.size()):
        feature = featuremap[i]
        chs = feature.getConvexHulls()
        for j in range(len(chs)):
            pts = chs[j].getHullPoints()
            hulls.loc[len(hulls)] = [pts.min(0)[0], pts.max(0)[0], pts.min(0)[1], pts.max(0)[1], False, -1]
    return hulls

def FeatureFindingMetabo1(mzfile):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(mzfile, exp)
    
    mtd_params = pyopenms.MassTraceDetection().getDefaults()
    mtd = pyopenms.MassTraceDetection()
    mtd.setParameters(mtd_params)
    mass_traces=[]
    mtd.run(exp, mass_traces)
    
    epdet_params = pyopenms.ElutionPeakDetection().getDefaults()
    epdet = pyopenms.ElutionPeakDetection()
    epdet.setParameters(epdet_params)
    splitted_mass_traces = []
    epdet.detectPeaks(mass_traces, splitted_mass_traces)
    
    ffm_params = pyopenms.FeatureFindingMetabo().getDefaults()
    ffm = pyopenms.FeatureFindingMetabo()
    ffm.setParameters(ffm_params)
    feature_map = pyopenms.FeatureMap()
    ffm.run(splitted_mass_traces, feature_map)
    return feature_map

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

def peaks2df(peaks):
    df = pd.DataFrame(columns=['rt', 'mz', 'intensity'])
    for i in range(peaks.shape[0]):
        rt  = peaks[i,3]
        mz  =peaks[i,0]
        intensity = peaks[i,6]
        df.loc[len(df)] = [rt, mz, intensity] 
    return df

def match_features(ground_truths, df):
    for i in range(len(df)):
        rt  = df.at[i, 'rt']
        mz  = df.at[i, 'mz']
        for j in range(len(ground_truths)):
            if(rt >= ground_truths.at[j, 'rt_min'] and rt <= ground_truths.at[j, 'rt_max'] and
               mz >= ground_truths.at[j, 'mz_min']-0.01 and mz <= ground_truths.at[j, 'mz_max']+0.01
               ):
                ground_truths.at[j, 'detected'] = True
                ground_truths.at[j, 'pic_id'] = i


if __name__=="__main__":
    mm48_all = pd.read_csv('simulation/MM48_annotations.csv')
    mm48_all['charge'] = [1] * mm48_all.shape[0]
    mm48_all['shape'] = ['gauss'] * mm48_all.shape[0]
    mm48_all['source'] = ['ESI'] * mm48_all.shape[0]
    mm48 = mm48_all[['Name', 'Formel','RT','RT2','Intensity','charge','shape','source']]
    mm48.to_csv('simulation/MM48_MSSimulator.csv', header=False, index=False)
    
    
    parameters = [[0, 1],[1, 10], [3,30], [5, 50], [10, 100],  [12, 120], [15, 150], [20, 200]]
    
    simulation('simulation/test.fasta','simulation/MM48_MSSimulator.csv',
               'MM48_MSS_Profile.mzML', 'MM48_MSS.featureXML', parameters[i][0]) 
    peak_picker = 'C:/Program Files/OpenMS/bin/PeakPickerHiRes.exe'
    subprocess.call([peak_picker,'-in', 'MM48_MSS_Profile.mzML',
                    '-out', 'MM48_MSS.mzML'])
    data2mzxml('MM48_MSS.mzML')
    ground_truths = parse_featureXML_GT('MM48_MSS.featureXML')
    
    
    
    mzfile =  "MM48_MSS.mzxml"
    mzMLfile =  "MM48_MSS.mzML"
    
    

    
    tic()
    feature_map = FeatureFindingMetabo(mzMLfile, parameters[i][1], 3)
    df_ffm = parse_featureXML_FFM(feature_map)
    toc()
    
    match_ffm = ground_truths.copy()
    match_features(match_ffm, df_ffm)
    match_ffm.detected.value_counts()
       
    from r_functions import XCMS
    tic()
    df_xcms = XCMS(mzMLfile, w1=5, w2=50, snr=3, intensity=parameters[i][1])
    toc()
    match_xcms = ground_truths.copy()
    match_features(match_xcms, df_xcms)
    match_xcms.detected.value_counts()
    
    tic()
    parser=mzXMLParser()
    lcms = parser.parseFile(mzfile.encode(sys.getfilesystemencoding()))
    pics_c = pm.FPICs(lcms, parameters[i][1], 200.0, 0.5)
    toc()
    match_fpic = ground_truths.copy()
    df_fpic = pics2df(pics_c)
    match_features(match_fpic, df_fpic)
    match_fpic.detected.value_counts()
    
    from FPIC import pics2peaks, merge_peaks
    df_fpic_merge = peaks2df(merge_peaks(pics2peaks(pics_c), 0.02, 10))
    match_fpic_merge = ground_truths.copy()
    match_features(match_fpic_merge, df_fpic_merge)
    match_fpic_merge.detected.value_counts()
    