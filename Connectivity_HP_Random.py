from brian2 import *
from scipy import sparse,linalg
from scipy.io import savemat
import random
seed(20)

from Parameters import N_PYh,N_INh, kPPh,kPIh,kIPh,kIIh,sigma_PYh,sigma_INh

G = N_PYh+N_INh     # total number of cells

C_PP_h = np.zeros((N_PYh,N_PYh))                     # append memory for connectivity matrices
C_PI_h = np.zeros((N_PYh,N_INh))
C_IP_h = np.zeros((N_INh,N_PYh))
C_II_h = np.zeros((N_INh,N_INh))

d_PP_h = np.zeros((N_PYh,N_PYh))
d_PI_h = np.zeros((N_PYh,N_INh))
d_IP_h = np.zeros((N_INh,N_PYh))
d_II_h = np.zeros((N_INh,N_INh))

cl_PP = np.zeros(N_PYh)
cl_PI = np.zeros(N_PYh)
cl_IP_P = np.zeros(N_INh)
cl_IP_I = np.zeros(N_INh)
cl_II = np.zeros(N_INh)

Ideg_E = abs(around(kIPh+(kIPh*0.05)*randn(N_INh)))  # E-degrees of all IN cells (abs in case of negative degree)
Edeg_E = abs(around(kPPh+(kPPh*0.05)*randn(N_PYh)))  # E-degrees of all PY cells 
Edeg_I = kPIh*ones(N_PYh)                            # I-degrees of all PY cells 
Ideg_I = abs(around(kIIh+(kIIh*0.05)*randn(N_INh)))  # I-degrees of all IN cells (abs in case of negative degree)

I = dot(11,range(1,int(G/11)+1))-1                # indexes of the IN cells
E = delete(range(G),I)                            # indexes of the PY cells
#-------------------------------------------------------------------------------

for i in range(N_PYh):                          # For each PY cell...
# E-E connections------------------------------
    CPP = delete(range(N_PYh),i)                # keep the numbers of all OTHER PY cells
    d_PP_h[i,:] = (E-E[i])                      # keep distances of E cells 
    index = random.sample(range(0, N_PYh-1), int(Edeg_E[i]))
    C_PP_h[i, CPP[index]] += 1
# E-I connections-------------------------------
    d_PI_h[i,:] = (I-E[i])                      # keep distances of I cells 
    index = random.sample(range(0, N_INh), int(Edeg_I[i]))
    C_PI_h[i, index] += 1
#-------------------------------------------------------------------------------
for i in range(N_INh):
# I-E connections-----------------------------
    d_IP_h[i,:] = abs(E-I[i])
    index = random.sample(range(0, N_PYh), int(Ideg_E[i]))
    C_IP_h[i, index] += 1
# I-I connections-------------------------------
    CII = delete(range(N_INh),i)
    d_II_h[i,:] = (I-I[i])                      # keep distances of I cells
    index = random.sample(range(0, N_INh-1), int(Ideg_I[i]))
    C_II_h[i, CII[index]] += 1
#-------------------------------------------------------------------------------
savemat('C_PP_h.mat', {'C_PP_h':array(C_PP_h)}) 
savemat('C_PI_h.mat', {'C_PI_h':array(C_PI_h)}) 
savemat('C_IP_h.mat', {'C_IP_h':array(C_IP_h)}) 
savemat('C_II_h.mat', {'C_II_h':array(C_II_h)}) 

savemat('d_PP_h.mat', {'d_PP_h':array(d_PP_h)}) 
savemat('d_PI_h.mat', {'d_PI_h':array(d_PI_h)}) 
savemat('d_IP_h.mat', {'d_IP_h':array(d_IP_h)}) 
savemat('d_II_h.mat', {'d_II_h':array(d_II_h)}) 