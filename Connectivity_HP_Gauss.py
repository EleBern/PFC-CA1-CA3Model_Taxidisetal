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
    #cl_PP[i] = np.count_nonzero(abs(d_PI_3[i,:])<=3*sigma_PY3)
    cl_PP[i] = np.sum(abs(d_PP_h[i,:])<=3*sigma_PYh)
    A = delete(d_PP_h[i,:],i)                   # remove the i cell itself
    A = exp(-A**2/(2*sigma_PYh**2))/(sqrt(2*pi*sigma_PYh**2)) # find the connection probability of each cell
    A = A/sum(A)                                # normalize probabilities (so they add up to 1)
    for j in range(1,A.size):                   # transform the probabilites into...
        A[j] = A[j] + A[j-1]                    # upper bounds of intervals of the (0,1) space
    conn = rand(int(Edeg_E[i]))                 # draw Edeg_E[i] random numbers from the [0,1) space
    for k in range(conn.size):                  # for each one (to find the interval it corresponds to)...
        z = A-conn[k]                           # subtract it from the upper bounds of the intervals
        z = np.asarray(np.where(z>=0))          # the first positive bound is the interval (and the cell)
        C_PP_h[i,CPP[z[0,0]]] += 1              # add a connetion to the right entry of the conn. matrix
      
# E-I connections-------------------------------
    d_PI_h[i,:] = (I-E[i])                      # keep distances of I cells 
    # cl_PI[i] = np.count_nonzero(abs(d_PI_3[i,:])<=3*sigma_PY3)
    cl_PI[i] = np.sum(abs(d_PI_h[i,:])<=3*sigma_PYh)
    A = exp(-d_PI_h[i,:]**2/(2*sigma_PYh**2))/(sqrt(2*pi*sigma_PYh**2)) # find the connection probabilities
    A = A/sum(A)                                # normalize
    for j in range(1,A.size):                   # transform them into intervals
        A[j] = A[j] + A[j-1]
    conn = rand(int(Edeg_I[i]))                 # draw Edeg_I[i] random numbers
    for k in range(conn.size):                  # find the interval they correspond to
        z = A-conn[k]
        z = np.asarray(np.where(z>=0)) 
        C_PI_h[i,z[0,0]] += 1                   # add connections to the appropriate entries.
#-------------------------------------------------------------------------------
for i in range(N_INh):
# I-E connections-----------------------------
    d_IP_h[i,:] = abs(E-I[i])
    # cl_IP_P[i] = np.count_nonzero(abs(d_PI_3[i,:])<=3*sigma_PY3)
    # cl_IP_I[i] = np.count_nonzero(abs(d_PI_3[i,:])<=3*sigma_IN3)
    cl_IP_P[i] = np.sum(abs(d_PI_h[i,:])<=3*sigma_PYh)
    cl_IP_I[i] = np.sum(abs(d_PI_h[i,:])<=3*sigma_INh)
    A = zeros(N_PYh)
    # near = find(abs(d_IP_3[i,:])<=sigma_IN3)
    near = np.nonzero(abs(d_IP_h[i,:])<=sigma_INh)
    A[near] = 1./np.count_nonzero(near) 
    A = A/sum(A)
    # Do not uncomment
#    A = exp(-d_IP_3[i,:]**2/(2*sigma_INh**2))/(sqrt(2*pi*sigma_INh**2))
#    A = A/sum(A)
    for j in range(1,A.size):
        A[j] = A[j] + A[j-1]
    conn = rand(int(Ideg_E[i]))
    for k in range(conn.size):
        z = A-conn[k]
        z = np.asarray(np.where(z>=0)) 
        C_IP_h[i,z[0,0]] += 1  
## I-I connections-------------------------------
    CII = delete(range(N_INh),i)
    d_II_h[i,:] = (I-I[i])                      # keep distances of I cells 
    #cl_II[i] = np.count_nonzero(abs(d_II_1[i,:])<=3*sigma_IN1)
    cl_II[i] = np.sum(abs(d_II_h[i,:])<=3*sigma_INh)
    A = delete(d_II_h[i,:],i)                   # remove the i cell itself
    A = exp(-A**2/(2*sigma_INh**2))/(sqrt(2*pi*sigma_INh**2))
    A = A/sum(A)
    for j in range(1,A.size):
        A[j] = A[j] + A[j-1]
    conn = rand(int(Ideg_I[i]))
    for k in range(conn.size):
        z = A-conn[k]
        z = np.asarray(np.where(z>=0))
        C_II_h[i,CII[z[0,0]]] += 1
#-------------------------------------------------------------------------------
savemat('C_PP_h.mat', {'C_PP_h':array(C_PP_h)}) 
savemat('C_PI_h.mat', {'C_PI_h':array(C_PI_h)}) 
savemat('C_IP_h.mat', {'C_IP_h':array(C_IP_h)}) 
savemat('C_II_h.mat', {'C_II_h':array(C_II_h)}) 

savemat('d_PP_h.mat', {'d_PP_h':array(d_PP_h)}) 
savemat('d_PI_h.mat', {'d_PI_h':array(d_PI_h)}) 
savemat('d_IP_h.mat', {'d_IP_h':array(d_IP_h)}) 
savemat('d_II_h.mat', {'d_II_h':array(d_II_h)}) 