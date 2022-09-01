from brian2 import *
from scipy import sparse,linalg
from scipy.io import savemat
import random
seed(20)

from Parameters import N_PYc,N_INc,N_PYh,N_INh,khcP,khcI,sigma_hcP,sigma_hcI

DoubleN_PYc = N_PYc * 2
DoubleN_INc = N_INc * 2
C_hc_P = zeros((N_PYh,N_PYc))                      
C_hc_I = zeros((N_PYh,N_INc))

deg_hcP = abs(around(khcP+(khcP*0.05)*randn(N_PYh)))      
deg_hcI = abs(around(khcI+(khcI*0.05)*randn(N_PYh)))

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
Gh = range(Gh)
#-------------------------------------------------------------------------------
  
for i in range(N_PYh):                          # For each PY cell in PFC...
# E-E connections------------------------------
    index = random.sample(range(0, N_PYc), int(deg_hcP[i]))
    C_hc_P[i, index] += 1
# E-I connections-------------------------------
    index = random.sample(range(0, N_INc), int(deg_hcI[i]))
    C_hc_I[i, index] += 1
#-------------------------------------------------------------------------------
savemat('C_hc_I.mat', {'C_hc_I':array(C_hc_I)}) 
savemat('C_hc_P.mat', {'C_hc_P':array(C_hc_P)}) 