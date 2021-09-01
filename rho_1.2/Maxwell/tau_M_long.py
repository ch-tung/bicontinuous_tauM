# -*- coding: utf-8 -*-
"""
Created on Sun May 16 10:31:20 2021
Calculate the maxwell relaxation time
@author: CHTUNG
"""

import numpy as np
from scipy.io import savemat
from scipy.io import loadmat
import time
tStart = time.time()

rho = 1.2
N = 10976

T = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0])
n_T = T.size
n_conf = 1
n_time = 1000000
n_correlation = 100000
n_inteval = 10

t = np.arange(0,n_correlation*n_inteval,n_inteval)
C_c = np.zeros((n_correlation,n_T));
G_V_c = np.zeros((n_T))
v = np.zeros((n_T))

for i_T in range(n_T):
    T_i = T[i_T]
    
    name_T = str('%.1f' % T_i)
    print(T_i)
    
    filename_thermo = './thermo' + name_T + '.mat'
    thermo_dict = loadmat(filename_thermo)
    thermo = thermo_dict['thermo']
      
    for i_c in range(n_conf):
        
        C = np.zeros((n_correlation))
        S = thermo[:,5:8,i_c]
        S_mean = np.average(S,0)
        S = S-(S_mean)
        
        G_V = np.average(np.average(S**2)/T_i)
        G_V_c[i_T] = G_V_c[i_T] + G_V
        
        for i_t in range(n_correlation):
            dt = i_t; # frames
            S_0 = S[0:n_time-dt,:]
            S_t = S[dt:n_time,:]
            C[i_t] = C[i_t] + np.average(np.average(S_0*S_t,0)/np.average(S**2,0))
            
        C_c[:,i_T] = C_c[:,i_T] + C
        
    C_c[:,i_T] = C_c[:,i_T]/n_conf
    G_V_c[i_T] = G_V_c[i_T]/n_conf
    
    v[i_T] = np.sqrt(G_V_c[i_T]/(rho*N))

    filename_maxwell = 'maxwell_long_py_' + name_T + '.mat'
    mdic = {'t':t, 'C_c':C_c, 'G_V_c':G_V_c, 'v':v}
    savemat(filename_maxwell, mdic)

tEnd = time.time()
print("It cost %f sec" % (tEnd - tStart))    
            