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


inf = gzip.open('F:/resources/isotope/pubchem/ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/2016-12-01/SDF/Compound_000000001_000025000.sdf.gz')
gzsuppl = Chem.ForwardSDMolSupplier(inf)
ms = [Chem.AddHs(x) for x in gzsuppl if x is not None]
ms_filtered = filter_pubchem(ms)
