from brian2 import *
from scipy import sparse,linalg
from scipy.io import savemat
import random
seed(20)

from Parameters import N_PYc,N_INc,N_PYh,N_INh,kchP,kchI,sigma_chP,sigma_chI

DoubleN_PYc = N_PYc * 2
C_ch_P = zeros((N_PYc,N_PYh))                      
C_ch_I = zeros((N_PYc,N_INh))

deg_chP = abs(around(kchP+(kchP*0.05)*randn(N_PYc)))      
deg_chI = abs(around(kchI+(kchI*0.05)*randn(N_PYc)))


Gc = N_PYc + N_INc
Ic = (dot(8,range(1,int(Gc/10)+1))-1)             # indexes of the IN cells
Ec = delete(range(Gc-len(Ic)),Ic)                        # indexes of the PY cells
Ec = Ec*2
Ic = (dot(9,range(1,int(Gc/5)+1))-1)         # indexes of the IN cells
for i in range(0,len(Ic),2):
    Ic[i] = Ic[i] - 1

Gh = N_PYh+N_INh
Ih = dot(11,range(1,int(Gh/11)+1))-1               # indexes of the IN cells
Eh = delete(range(Gh),Ih)                          # indexes of the PY cells
#-------------------------------------------------------------------------------
  
for i in range(N_PYc):                          # For each PY cell in PFC...
# E-E connections------------------------------
    index = random.sample(range(0, N_PYh), int(deg_chP[i]))
    C_ch_P[i, index] += 1
# E-I connections-------------------------------
    index = random.sample(range(0, N_INh), int(deg_chI[i]))
    C_ch_I[i, index] += 1

#-------------------------------------------------------------------------------
savemat('C_ch_P.mat', {'C_ch_P':array(C_ch_P)}) 
savemat('C_ch_I.mat', {'C_ch_I':array(C_ch_I)}) 