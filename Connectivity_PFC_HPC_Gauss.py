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
    d_chP = Eh-Ec[i]   # WORKS ONLY IF N_PYc = N_PYh!!!!!!    # keep distances of all HPC PY cells     
    A = exp(-d_chP**2/(2*sigma_chP**2))/(sqrt(2*pi*sigma_chP**2)) # find the connection probability of each cell
    A = A/sum(A)                                # normalize probabilities (so they add up to 1) 
    for j in range(1,A.size):                   # transform the probabilites into...
        A[j] = A[j] + A[j-1]                    # upper bounds of intervals of the (0,1) space
    conn = rand(int(deg_chP[i]))                # draw deg_11[i] random numbers from the [0,1) space
    for k in range(conn.size):                  # for each one (to find the interval it corresponds to)...
        z = A-conn[k]                           # subtract it from the upper bounds of the intervals
        z = np.asarray(np.where(z>=0))          # the first positive bound is the interval (and the cell)
        C_ch_P[i,z[0,0]] += 1 

# E-I connections------------------------------
    d_chI = Ih-Ec[i]   # WORKS ONLY IF N_PYc = N_PYh!!!!!!    # keep distances of all HPC IN cells     
    A = exp(-d_chI**2/(2*sigma_chI**2))/(sqrt(2*pi*sigma_chI**2)) # find the connection probability of each cell
    A = A*1./sum(A)                             # normalize probabilities (so they add up to 1) 
    for j in range(1,A.size):                   # transform the probabilites into...
        A[j] = A[j] + A[j-1]                    # upper bounds of intervals of the (0,1) space
    conn = rand(int(deg_chI[i]))                # draw deg_11[i] random numbers from the [0,1) space
    for k in range(conn.size):                  # for each one (to find the interval it corresponds to)...
        z = A-conn[k]                           # subtract it from the upper bounds of the intervals
        z = np.asarray(np.where(z>=0))          # the first positive bound is the interval (and the cell)
        C_ch_I[i,z[0,0]] += 1
#-------------------------------------------------------------------------------
savemat('C_ch_P.mat', {'C_ch_P':array(C_ch_P)}) 
savemat('C_ch_I.mat', {'C_ch_I':array(C_ch_I)}) 