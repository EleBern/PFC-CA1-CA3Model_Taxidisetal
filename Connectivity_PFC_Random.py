from brian2 import *
from scipy import sparse,linalg
from scipy.io import savemat
import random
seed(20)

from Parameters import N_PYc,N_INc,Gc, kPPc,kPIc,kIPc,kIIc, sigma_PYc,sigma_INc

CPP = np.zeros((N_PYc,N_PYc))                       # append memory for connectivity matrices
CPI = np.zeros((N_PYc,N_INc))
CIP = np.zeros((N_INc,N_PYc))
CII = np.zeros((N_INc,N_INc))

Ideg_E = abs(around(kIPc+5.*randn(N_INc)))          # E-degrees of all IN cells (abs in case of negative degree)
Ideg_I = abs(around(kIIc+5.*randn(N_INc)))          # I-degrees of all IN cells 
Edeg_E = abs(around(kPPc+5.*randn(N_PYc)))          # E-degrees of all PY cells 
Edeg_I = abs(around(kPIc+5.*randn(N_PYc)))          # I-degrees of all PY cells 

I = dot(5,range(1,int(Gc/5)+1))-1                   # indexes of the IN cells
E = delete(range(Gc),I)                             # indexes of the PY cells
#-------------------------------------------------------------------------------

for i in range(N_PYc):                          # For each PY cell...
# E-E connections------------------------------
    CP = delete(range(N_PYc),i)                 # keep the numbers of all OTHER PY cells
    prova = random.sample(range(0, N_PYc-1), int(Edeg_E[i]))
    CPP[i, CP[prova]] += 1
# E-I connections-------------------------------
    prova = random.sample(range(0, N_INc), int(Edeg_I[i]))
    CPI[i, prova] += 1
#-------------------------------------------------------------------------------

for i in range(N_INc):
# I-I connections-------------------------------
    CI = delete(range(N_INc),i)
    prova = random.sample(range(0, N_INc-1), int(Ideg_I[i]))
    CII[i, CI[prova]] += 1
# I-E connections------------------------------
    prova = random.sample(range(0, N_PYc), int(Ideg_E[i]))
    CIP[i, prova] += 1
#-------------------------------------------------------------------------------
savemat('CIP.mat', {'CIP':array(CIP)}) 
savemat('CII.mat', {'CII':array(CII)}) 
savemat('CPP.mat', {'CPP':array(CPP)}) 
savemat('CPI.mat', {'CPI':array(CPI)}) 
#-------------------------------------------------------------------------------