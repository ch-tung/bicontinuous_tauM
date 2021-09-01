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
dr = 0.25
rmax = 10
rr = np.linspace(dr, rmax, int(rmax/dr), endpoint=True)
len_rr = rr.size

dq = 0.1*np.pi
qmax = 4*np.pi
qq = np.linspace(dq, qmax, int(qmax/dq), endpoint=True)
len_qq = qq.size

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
    
    return r_all, l, v_all

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

###############################################
#%%
#Temperature = np.arange(0,20,1)*0.1 + 0.1
#Temperature = np.arange(0,3,1)+3
Temperature = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0])
#Temperature = np.array([1.0])

n_T = Temperature.shape[0]
n_t = t.shape[0]

for iT in range(n_T):
    T_i = Temperature[iT]
    name_T = str('%.1f' % T_i)
    print(T_i)
        
    sim_r = np.zeros((len_rr,n_frame))
    sim_q = np.zeros((len_qq,n_frame))
    sim_s = np.zeros((1,n_frame))
    sim_a = np.zeros((1,n_frame))
    sim_r_c = np.zeros((len_rr,n_frame))
    sim_q_c = np.zeros((len_qq,n_frame))
    sim_s_c = np.zeros((1,n_frame))
    sim_a_c = np.zeros((1,n_frame))
        
    for ic in range(n_conf):
        name_P = str(ic+1)
        print(ic+1)
        
        filename_dump = '../' + name_T + '/eq.' + name_P + '.dump'
        r_all, l, v_all = loaddump(filename_dump)
        L = l[:,1]-l[:,0]
        
        # time origin
        it = 0
        v0 = v_all[:,:,it]
        p0 = np.sum(v0[:,0:3],1)/3
        s0 = v0[:,3];
        
        v_all_t = v_all[:,:,it]
        W_0, V_0 = calculate_shape(v_all_t,n_particle)
        V_0 = np.transpose(V_0,[0, 2, 1])
                
        rt = r_all[:,:,it]
        sl = np.arange(0,n_particle,1)
        index_bin = dist_index_matrix(rt, L, rc, sl)
        index_bin_triu = np.triu(index_bin)
        
        for it in range(n_t):
            vt = v_all[:,:,it]
            pt = np.sum(vt[:,0:3],1)/3
            st = vt[:,3];
            
            v_all_t = v_all[:,:,it]
            W, V = calculate_shape(v_all_t,n_particle)
            V = np.transpose(V,[0, 2, 1])
            ###############################################
            #%% calculate correlation
            Cv = np.zeros((n_particle,3))
            
            for iV in range(3):
                V_0_i = V_0[:,:,iV];
                V_i = V[:,:,iV];
                Cv[:,iV] = np.sum(np.multiply(V_0_i,V_i),0)**2
                
            Cv = (Cv*3-1)/2
            Cvz = (Cv-0*np.average(Cv,0))
            
            ###############################################
            #%% calculate similiarity
            iV = 2
            sim = Cvz[:,iV].reshape(n_particle,1)*Cvz[:,iV].reshape(1,n_particle)
            sim_triu = np.triu(sim)
            
            for i_index in range(len_rr):
                sim_r[i_index,it] = np.average(sim[index_bin==i_index+1])
                
            #i_NaN = np.isnan(sim_r[:,it])
            #sim_r[i_NaN,it] = 0
            #sim_r[rr<1,it] = 0
            
            for i_index in range(len_qq):
                sim_q[i_index,it] = np.sum(sim_r[:,it]*rr*np.sin(rr*qq[i_index]))/qq[i_index]
            
            # self
            sim_s[:,it] = np.average(np.diag(sim))
            # all
            sim_a[:,it] = np.average(Cvz[:,iV],0)**2
            
        sim_r_c = sim_r_c + sim_r
        sim_q_c = sim_q_c + sim_q
        sim_s_c = sim_s_c + sim_s
        sim_a_c = sim_a_c + sim_a
        
    sim_r_c = sim_r_c/n_conf
    sim_q_c = sim_q_c/n_conf
    sim_s_c = sim_s_c/n_conf
    sim_a_c = sim_a_c/n_conf
    
    filename_dcorr = 'qcorr4_vo_py_' + name_T + '.mat'
    mdic = {'sim_r_c':sim_r_c, 'sim_q_c':sim_q_c, 'sim_s_c':sim_s_c, 'sim_a_c':sim_a_c, 'rr':rr, 'qq':qq, 't':t}
    savemat(filename_dcorr, mdic)
           
###############################################
#%% calculate similiarity
tEnd = time.time()
print("It cost %f sec" % (tEnd - tStart))