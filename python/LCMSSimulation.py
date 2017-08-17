# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 08:14:15 2017

@author: zmzhang
"""


import subprocess
import pyopenms


def simulation(fasta, contaminants, out, out_cntm,
               simulator = 'C:/Program Files/OpenMS/bin/MSSimulator.exe'):   
    """
        Should copy "C:\Program Files\OpenMS\share\OpenMS\examples" to working directory of Python
    """ 
   
    subprocess.call([simulator, '-in', fasta, '-out', out, '-out_cntm',out_cntm, 
               '-algorithm:MSSim:RawSignal:contaminants:file', contaminants,
               '-algorithm:MSSim:Ionization:mz:lower_measurement_limit', '10',
               '-algorithm:MSSim:Ionization:mz:upper_measurement_limit', '1000',
               '-algorithm:MSSim:RT:total_gradient_time', '1000',
               '-algorithm:MSSim:RT:sampling_rate', '0.1',
               '-algorithm:MSSim:RT:scan_window:min', '0',
               '-algorithm:MSSim:RT:scan_window:max', '1000'])

def parse_featureXML(feature_file):
    featuremap = pyopenms.FeatureMap()
    featurexml = pyopenms.FeatureXMLFile()
    featurexml.load(feature_file, featuremap)
    
    for i in range(featuremap.size()):
        feature = featuremap[i]
        chs = feature.getConvexHulls()
        pts = chs[0].getHullPoints()
    return pts

if __name__=="__main__":
    import pandas as pd
    mm48_all = pd.read_csv('simulation/MM48_annotations.csv')
    mm48_all['charge'] = [1] * mm48_all.shape[0]
    mm48_all['shape'] = ['gauss'] * mm48_all.shape[0]
    mm48_all['source'] = ['ESI'] * mm48_all.shape[0]
    mm48 = mm48_all[['Name', 'Formel','RT','RT2','Intensity','charge','shape','source']]
    mm48.to_csv('simulation/MM48_MSSimulator.csv', header=False, index=False)
    
    simulation('simulation/test.fasta','simulation/MM48_MSSimulator.csv',
               'MM48_MSS.mzML', 'MM48_MSS.featureXML' )
    
    parse_featureXML('MM48_MSS.featureXML')
    
