from brian2 import *

# mV, cm, ms, msiemens/cm**2, uA/cm**2, ufarad/cm**2, mmole   

# NETWORK & CELL SIZES (in cm**2)-----------------------------------------------
N_PYc = 1000         
N_INc = 250
          
N_PY3 = 1000
N_IN3 = 100

N_PY1 = 1000
N_IN1 = 100

N_PYh = N_PY3 + N_PY1
N_INh = N_IN3 + N_IN1
N_IN  = N_INc + N_INh

Gc = N_PYc + N_INc     # total number of cells
L = 0.5                # total length of PFC network (5*mm)

As = 1.5e-4     # Dendrtic compartment size (0.015*(mm**2))
Ad = 3.5e-4     # Axosomatic compartment size (0.035*(mm**2)) 
S_PY = As + Ad  # Total pyramidal cell size (0.05*(mm**2))
S_IN = 2e-4     # Total interneuron size (0.02*(mm**2))
#-------------------------------------------------------------------------------

# CORTICAL PYRAMIDAL CELLS Parameters-------------------------------------------
C_c = 1.   
phi = 4.      
gNa_c = 50. 
gK_c = 10.5  
gA_c = 1.   
gKS_c = 0.576  
gKNa_c = 1.33      
gCa_c = 0.43    
gKCa_c = 0.57   
gNaP_c = 0.0686 
gAR_c = 0.0257  

alpha_Na = 10.    # 0.01*mmole*(nA**-1)*(ms**-1)  turned to mmole*(uA**-1)*(ms**-1)
R_pump = 0.018    # mmole*ms**-1
Na_eq = 9.5       # mmole
alpha_Ca = 0.005  # umole*(nA**-1)*(ms**-1)   same as  mmole*(uA**-1)*(ms**-1)
K_D = 0.03        # mmole
tau_hA = 15.
tau_Ca = 150. 

VL_c = -60.95 
gL_c = 0.0667
gsd_c = 1.75e-3   # msiemens

VK_c  = -100.
VNa_c = 55.
VCa_c = 120.
#-------------------------------------------------------------------------------

# HIPPOCAMPAL PYRAMIDAL CELLS Parameters----------------------------------------
C_h = 3. 
p = 0.5
gNa_h = 30.
gK_h = 15. 
gKAHP_h = 0.8   
gKCa_h = 15.

gCa_3 = 10.
gCa_1 = 7.

VL_h = -60.
gL_h = 0.1
gsd_h = 2.1

VK_h  = -75.
VNa_h = 60.
VCa_h = 80.
# ------------------------------------------------------------------------------

# INTERNEURONS Parameters-------------------------------------------------------
Ci = 1.  
gNai = 35.
gKi = 9.     
VKi = -90.
VNai = 55.

VLi_m = -63.8 
gLi_m = 0.1025
# ------------------------------------------------------------------------------

# SYNAPTIC Parameters-----------------------------------------------------------
gAMPA_PP_c = 30* 5.4*(1.e-6/Ad)    #29.3, 29.5, 30
gAMPA_PI_c = 10*2.25*(1.e-6/S_IN)
gGABA_IP_c = 4.15*(1.e-6/As)  #1-6
gGABA_II_c = 10*0.165*(1.e-6/S_IN) 

gAMPA_PP_h = 50.*(1.e-6/S_PY)  
gAMPA_PI_h = 80.*(1.e-6/S_IN)  
gGABA_IP_h = 100.*(1.e-6/S_PY) 
gGABA_II_h = 40.*(1.e-6/S_IN)

# Destexhe synaptic parameters
# Taken from modelDB
# Thalamocortical and Thalamic Reticular Network (Destexhe et al 1996)
Cmax_GABA = 0.5 
Cdur_GABA = 0.3
Alpha_GABA = 10.5
Beta_GABA = 0.166
Erev_GABA = -80.
Deadtime_GABA = 1.

# AMPA
Cmax_AMPA = 0.5 
Cdur_AMPA = 0.3
# Alpha_AMPA = 0.94
Alpha_AMPA_h = 3.94
Alpha_AMPA_c = 1.54
Beta_AMPA = 0.18
Erev_AMPA = 0.
Deadtime_AMPA = 1.
# ------------------------------------------------------------------------------

# HPC AXONAL DELAYS---------------------------------------------------
dx = 10e-3         # 10 microns = 10*10^-3 mm
V_ax_P = 0.5       # PY axonal conductance velocity (mm/msec)
V_ax_I = 0.1       # IN axonal conductance velocity
#-------------------------------------------------------------------------------

# PFC CONNECTIVITY Parameters---------------------------------------------------
kPPc = kPIc = kIPc = kIIc = 20. 
sigma_PYc = 0.0250/(L/Gc)  # PY connectivity deviation in cell numbers
sigma_INc = 0.0125/(L/Gc)  # IN connectivity deviation in cell numbers
#-------------------------------------------------------------------------------

# HPC CONNECTIVITY Parameters---------------------------------------------------
kPPh = 200. 
kPIh = 20.   
kIPh = 400. 
kIIh = 100.
sigma_PYh = 100.   # PY connectivity deviation in cell numbers
sigma_INh = 30.    # IN connectivity deviation in cell numbers
#-------------------------------------------------------------------------------

# PFC-to-HPC CONNECTIVITY Parameters--------------------------------------------
# alpha_PFC_CAh_P = 2.#50.
# alpha_PFC_CAh_I = 2.
kchP = 50.                  
kchI = 2.
sigma_chP = 50.
sigma_chI = 50.
#-------------------------------------------------------------------------------

# HPC-to-PFC CONNECTIVITY Parameters--------------------------------------------
# alpha_CA1_PFC_P = 8.   
# alpha_CA1_PFC_I = 13.   
khcP = 50.
khcI = 2.        
sigma_hcP = 50.
sigma_hcI = 50.
#-------------------------------------------------------------------------------