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
    A = (E-E[i])                                # keep distances of E cells 
    A = delete(A,i)                             # remove the i cell itself
    A = exp(-A**2/(2*sigma_PYc**2))/(sqrt(2*pi*sigma_PYc**2)) # find the connection probability of each cell
    A = A/sum(A)                                # normalize probabilities (so they add up to 1)
    for j in range(1,A.size):                   # transform the probabilites into...
        A[j] = A[j] + A[j-1]                    # upper bounds of intervals of the (0,1) space
    conn = rand(int(Edeg_E[i]))                 # draw Edeg_E[i] random numbers from the [0,1) space
    for k in range(conn.size):                  # for each one (to find the interval it corresponds to)...
        z = A-conn[k]                           # subtract it from the upper bounds of the intervals
        z = np.asarray(np.where(z>=0))          # the first positive bound is the interval (and the cell)
        CPP[i,CP[z[0,0]]] += 1                  # add a conncetion to the right entry of the conn. nmatrix
# E-I connections-------------------------------
    A = (I-E[i])                                # keep distances of I cells 
    A = exp(-A**2/(2*sigma_PYc**2))/(sqrt(2*pi*sigma_PYc**2)) # find the connection probabilities
    A = A/sum(A)                                # normalize
    for j in range(1,A.size):                   # transform them into intervals
        A[j] = A[j] + A[j-1]
    conn = rand(int(Edeg_I[i]))                 # draw Edeg_I[i] random numbers
    for k in range(conn.size):                  # find the interval they correspond to
        z = A-conn[k]
        z = np.asarray(np.where(z>=0))
        CPI[i,z[0,0]] += 1                      # add connections to the appropriate entires.
#-------------------------------------------------------------------------------

for i in range(N_INc):
# I-I connections-------------------------------
    CI = delete(range(N_INc),i)
    A = (I-I[i])
    A = delete(A,i)
    A = exp(-A**2/(2*sigma_INc**2))/(sqrt(2*pi*sigma_INc**2))
    A = A/sum(A)
    for j in range(1,A.size):
        A[j] = A[j] + A[j-1]
    conn = rand(int(Ideg_I[i]))
    for k in range(conn.size):
        z = A-conn[k]
        z = np.asarray(np.where(z>=0))
        CII[i,CI[z[0,0]]] += 1
# I-E connections------------------------------
    A = (E-I[i])
    A = exp(-A**2/(2*sigma_INc**2))/(sqrt(2*pi*sigma_INc**2))
    A = A/sum(A)
    for j in range(1,A.size):
        A[j] = A[j] + A[j-1]
    conn = rand(int(Ideg_E[i]))
    for k in range(conn.size):
        z = A-conn[k]
        z = np.asarray(np.where(z>=0))
        CIP[i,z[0,0]] += 1
#-------------------------------------------------------------------------------
savemat('CIP.mat', {'CIP':array(CIP)}) 
savemat('CII.mat', {'CII':array(CII)}) 
savemat('CPP.mat', {'CPP':array(CPP)}) 
savemat('CPI.mat', {'CPI':array(CPI)}) 
