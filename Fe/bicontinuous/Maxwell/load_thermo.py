# -*- coding: utf-8 -*-
"""
Created on Mon May 10 08:28:47 2021
load thermo data from LAMMPS .out file
@author: CHTUNG
"""

import numpy as np
from scipy.io import savemat

filename_out = 'maxwell_long.out'
T = np.array([1500, 2300, 5000])

n_T = T.size
n_conf = 1
n_time = 1000001

header = np.zeros([n_conf*n_T])

with open(filename_out) as out, open('linenumber.txt', 'w') as line_out:
    keyword = 'Step Temp Pxx Pyy Pzz Pxy Pxz Pyz'
    linenumber = 0
    indexnumber = 0
    for line in out:
        linenumber = linenumber +1
        if keyword in line:
            header[indexnumber] = linenumber
            indexnumber = indexnumber +1
            print(linenumber, file=line_out)
            
with open(filename_out,'r') as fp:
    lines = fp.readlines()

for iT in range(n_T):
    T_i = T[iT]
    name_T = str('%d' % T_i)
    
    data = np.zeros([n_time,8,n_conf])
    for i_conf in range(n_conf):
        index_header = iT+n_T*i_conf
        print(index_header)
        hi = np.int(header[index_header])
        data[:,:,i_conf] = np.genfromtxt(lines[hi:hi+n_time])
        
    filename_thermo = 'thermo' + name_T + '.mat' #
    mdic = {'thermo':data}
    savemat(filename_thermo, mdic)
