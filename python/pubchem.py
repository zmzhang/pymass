# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 09:55:42 2017

@author: zmzhang
"""

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from rdkit.Chem.rdmolops import GetFormalCharge
import gzip
from collections import Counter
import json, glob, os

def filter_pubchem(ms):
    ms_filtered = []
    elements = set(['C', 'H', 'O', 'N', 'S', 'P', 'Cl', 'B','Br','Se'])
    for m in ms:
        mw = CalcExactMolWt(m)
        if mw < 100 or mw> 1500:
            continue
        
        if GetFormalCharge(m) != 0:
            continue
        
        atoms = [a.GetSymbol() for a in m.GetAtoms()]
        c = Counter(atoms)
        if 'C' in c and 'H' in c:
            if 'S' in c and c['S'] > 5:
                continue
            if 'Cl' in c and c['Cl'] > 5:
                continue 
            if 'Br' in c and c['Br'] > 5:
                continue
            if 'B' in c and c['B'] > 5:
                continue  
            if set(c.keys()).issubset(elements):
                ms_filtered.append(CalcMolFormula(m))
    return ms_filtered

def write_json(ms, fname):
    f = open(fname, 'w')
    json.dump(ms, f)
    f.close()


sdf_path = 'F:/resources/isotope/pubchem/ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/2016-12-01/SDF/'
sdfs =  glob.glob(sdf_path + '*.sdf.gz')
for i, sdf in enumerate(sdfs):
    print( '{:2.4f}%  {}'.format((100 * float(i)/len(sdfs)), sdf))
    if not(os.path.exists(sdf[0:-6]+'json')):
        sdf_gz = gzip.open(sdf)
        gzsuppl = Chem.ForwardSDMolSupplier(sdf_gz)
        ms = [Chem.AddHs(x) for x in gzsuppl if x is not None]
        ms_filtered = filter_pubchem(ms)
        write_json(ms_filtered, sdf[0:-6]+'json')


