"""
Mar 2021
calculate dynamical correlation length
@author: CHTUNG
"""
import numpy as np
import numpy.matlib
from scipy.io import loadmat
from scipy.io import savemat

import time
tStart = time.time()

###############################################
# define time interval
t = np.arange(1,10,1)
for it_p in range(7):
    t = np.append(t,np.arange(1,10,1)*(10**(it_p+1)))

t = t[np.where(t<=10000000)]
t = np.append(0,t)
t = t

###############################################
# define simulation
n_conf = 4
n_particle = 6912
n_frame = t.size
rs = 3
rc = 1
dr = 0.1
rmax = 10
rr = np.linspace(dr, rmax, int(rmax/dr), endpoint=True)
len_rr = rr.size

###############################################
#%%
# define functions
def is_header(x):
    particlenumber = n_particle
    return np.remainder(x,particlenumber+9)>8

def loaddump(filename):
    # load file
    with open(filename,'r') as fp:
        lines = fp.readlines()
    lines_header = lines[0:9]
    l = np.genfromtxt(lines_header[5:8], delimiter=' ')
    
    # remove headers
    particlenumber = n_particle
    nframe = n_frame
    index_all = range((particlenumber+9)*nframe)
    
    index = list(filter(is_header, index_all))
    
    from operator import itemgetter
    lines_mod = list(itemgetter(*index)(lines))
    
    # convert to array
    data_np = np.genfromtxt(lines_mod, delimiter=' ')
    
    r = data_np[:,2:5]
    
    # reshape and permute
    r_all = np.reshape(r,(nframe,particlenumber,3))
    r_all = r_all.transpose((1,2,0))
    
    v = data_np[:,5:11]
    
    # reshape and permute
    v_all = np.reshape(v,(nframe,particlenumber,6))
    v_all = v_all.transpose((1,2,0))

    Voro = data_np[:,14]
    
    return r_all, l, v_all, Voro

def dist_index_matrix(rt, L, rc, sl):
    Neig = sl
    N = Neig.size
    rjk = np.zeros((N,N,3))
    rjk = rt[sl,:].reshape(N,1,3) - rt[sl,:].reshape(1,N,3)
    
    for i_d in range(3):
        rjk[:,:,i_d] = rjk[:,:,i_d] - np.around(rjk[:,:,i_d]/L[i_d],0)*L[i_d]
    
    djk2 = np.sum(np.power(rjk,2),2)
    djk = np.sqrt(djk2)
    
    index_bin = np.ceil(djk/dr)
    
    return index_bin

def calculate_shape(v_all_t,n_particle):
    W = np.zeros((3,n_particle))
    V = np.zeros((3,3,n_particle))
    
    for i in range(n_particle):
        si = np.array([[v_all_t[i,0], v_all_t[i,3], v_all_t[i,4]],[v_all_t[i,3], v_all_t[i,1], v_all_t[i,5]],[v_all_t[i,4], v_all_t[i,5], v_all_t[i,2]]])
        Wi,Vi = np.linalg.eig(si)
        index_eig = np.argsort(Wi)
        W[:,i] = Wi[index_eig]
        V[:,:,i] = Vi[:,index_eig]       
        
    return W, V
	
def calculate_shape_rg(rg_all,n_particle):
    W = np.zeros((3,n_particle))
    V = np.zeros((3,3,n_particle))
    
    for i in range(n_particle):
        Wi,Vi = np.linalg.eig(rg_all[:,:,i])
        index_eig = np.argsort(Wi)
        W[:,i] = Wi[index_eig]
        V[:,:,i] = Vi[:,index_eig]       
        
    return W, V

def CG(rt, L, rc, sl, p):
    Neig = sl
    N = Neig.size
    rjk = np.zeros((N,N,3))
    rjk = rt[sl,:].reshape(N,1,3) - rt[sl,:].reshape(1,N,3)
    
    for i_d in range(3):
        rjk[:,:,i_d] = rjk[:,:,i_d] - np.around(rjk[:,:,i_d]/L[i_d],0)*L[i_d]
    
    djk2 = np.sum(np.power(rjk,2),2)
    
    weight = np.exp(-djk2/2/rc)
    weight_sum = np.sum(weight,1)
    
    p_CG = np.zeros(p.shape)
    
    for i_1 in range(N):
        p_CG[i_1,:] = np.sum((p.T*weight[:,i_1]).T,0)/weight_sum[i_1]
        
    return p_CG


###############################################
#%%
#Temperature = np.arange(0,20,1)*0.1 + 0.1
#Temperature = np.arange(0,3,1)+3
Temperature = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0])

n_T = Temperature.shape[0]
n_t = t.shape[0]

for iT in range(n_T):
    T_i = Temperature[iT]
    name_T = str('%.1f' % T_i)
    print(T_i)
        
    filename_rg = '../rg_mat/rg_py_' + name_T + '.mat'
    rg_all_dict = loadmat(filename_rg)
    rg_all_t_c = rg_all_dict['rg_all_t_c']
		
    C_0 = np.zeros(n_frame)
    C_1 = np.zeros(n_frame)
    
    C_0_c = np.zeros(n_frame)
    C_1_c = np.zeros(n_frame)
    
    Var_s = np.zeros(n_frame)
        
    for ic in range(n_conf):
        name_P = str(ic+1)
        print(ic+1)
        
        filename_dump = '../' + name_T + '/eq.' + name_P + '.dump'
        r_all, l, v_all, Voro = loaddump(filename_dump)
        L = l[:,1]-l[:,0]
        
        # time origin
        it = 0
        v0 = v_all[:,:,it]
        p0 = np.sum(v0[:,0:3],1)/3
        s0 = v0[:,3:6]
        Ms0 = np.average(s0,0)
        
        v_all_t = v_all[:,:,it]
        W_0, V_0 = calculate_shape(v_all_t,n_particle)
        V_0 = np.transpose(V_0,[0, 2, 1])
		
        rg_all = rg_all_t_c[:,:,:,it,ic]
        W_0_rg, V_0_rg = calculate_shape_rg(rg_all,n_particle)
        V_0_rg = np.transpose(V_0_rg,[0, 2, 1])  # structure eigenvector
                
        rt = r_all[:,:,it]
        sl = np.arange(0,n_particle,1)

        v0_CG = CG(rt, L, rc, sl, v0)
        p0_CG = np.sum(v0_CG[:,0:3],1)/3
        s0_CG = v0_CG[:,3];
        Ms0_CG = np.average(s0_CG,0)
        
        for it in range(n_t):
            vt = v_all[:,:,it]
            pt = np.sum(vt[:,0:3],1)/3
            st = vt[:,3:6]
            Mst = np.average(st,0)
            
            v_all_t = v_all[:,:,it]
            W, V = calculate_shape(v_all_t,n_particle)
            V = np.transpose(V,[0, 2, 1])
			
            rg_all = rg_all_t_c[:,:,:,it,ic]
            W_rg, V_rg = calculate_shape_rg(rg_all,n_particle)
            V_rg = np.transpose(V_rg,[0, 2, 1])  # structure eigenvector
            
            rt = r_all[:,:,it]
            
            vt_CG = CG(rt, L, rc, sl, vt)
            pt_CG = np.sum(vt_CG[:,0:3],1)/3
            st_CG = vt_CG[:,3];
            Mst_CG = np.average(st_CG,0)
			
            ###############################################
            #%% calculate correlation
            Cs = ((s0_CG)*(st_CG))
			
            Cv_rg = np.zeros((n_particle,3))
            
            for iV in range(3):
                V_0_rg_i = V_0_rg[:,:,iV];
                V_rg_i = V_rg[:,:,iV];
                Cv_rg[:,iV] = np.sum(np.multiply(V_0_rg_i,V_rg_i),0)**2
                
            Cv_rg = (Cv_rg*3-1)/2
            Cvz_rg = (Cv_rg-np.average(Cv_rg,0))
            
            C_0[it] = np.average(Cs[Cvz_rg[:,2]>=0])
            C_1[it] = np.average(Cs[Cvz_rg[:,2]<0])
            
        C_0_c = C_0_c + C_0
        C_1_c = C_1_c + C_1
        Var_s = Var_s + np.average((s0_CG)**2)
        
    C_0_c = C_0_c/n_conf
    C_1_c = C_1_c/n_conf
    Var_s = Var_s/n_conf
    
    filename_dcorr = 'tcorr_CG_vs_py_' + name_T + '.mat'
    mdic = {'C_0':C_0_c, 'C_1':C_1_c, 't':t, 'Var_s':Var_s}
    savemat(filename_dcorr, mdic)
           
###############################################
#%% calculate similiarity
tEnd = time.time()
print("It cost %f sec" % (tEnd - tStart))