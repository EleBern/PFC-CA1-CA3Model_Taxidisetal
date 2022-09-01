from brian2 import *

# mV, cm, ms, msiemens/cm**2, uA/cm**2, ufarad/cm**2, mmole   

# NETWORK & CELL SIZES (in cm**2)-----------------------------------------------          
N_PY3 = 1000
N_IN3 = 100

N_PY1 = 1000
N_IN1 = 100

N_PYh = N_PY3 + N_PY1
N_INh = N_IN3 + N_IN1
N_IN  = N_INh


As = 1.5e-4     # Dendrtic compartment size (0.015*(mm**2))
Ad = 3.5e-4     # Axosomatic compartment size (0.035*(mm**2)) 
S_PY = As + Ad  # Total pyramidal cell size (0.05*(mm**2))
S_IN = 2e-4     # Total interneuron size (0.02*(mm**2))
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
# gAMPA_PP_3 = 1.*(1.e-6/S_PY)  
# gAMPA_PI_3 = 1.*(1.e-6/S_IN)  
# gGABA_IP_3 = 1.*(1.e-6/S_PY) 
# gAMPA_PP_1 = 1.*(1.e-6/S_PY) 
# gAMPA_PI_1 = 1.*(1.e-6/S_IN)  
# gGABA_IP_1 = 1.*(1.e-6/S_PY)
# gGABA_II_1 = 1.*(1.e-6/S_IN)

# gAMPA_PP_3 = 70.*(1.e-6/S_PY)  
# gAMPA_PI_3 = 70.*(1.e-6/S_IN)  
# gGABA_IP_3 = 70.*(1.e-6/S_PY) 
# gAMPA_PP_1 = 70.*(1.e-6/S_PY) 
# gAMPA_PI_1 = 70.*(1.e-6/S_IN)  
# gGABA_IP_1 = 70.*(1.e-6/S_PY)
# gGABA_II_1 = 70.*(1.e-6/S_IN)

# gAMPA_PP_3 = 50.*(1.e-6/S_PY)  
# gAMPA_PI_3 = 80.*(1.e-6/S_IN)  
# gGABA_IP_3 = 100.*(1.e-6/S_PY) 
# gAMPA_PP_1 = 80.*(1.e-6/S_PY) 
# gAMPA_PI_1 = 80.*(1.e-6/S_IN)  
# gGABA_IP_1 = 80.*(1.e-6/S_PY)
# gGABA_II_1 = 40.*(1.e-6/S_IN)

gAMPA_PP_h = 50.*(1.e-6/S_PY)  
gAMPA_PI_h = 80.*(1.e-6/S_IN)  
gGABA_IP_h = 100.*(1.e-6/S_PY) 
gGABA_II_h = 40.*(1.e-6/S_IN)

# Destexhe synaptic parameters
# Taken from modelDB
# Thalamocortical and Thalamic Reticular Network (Destexhe et al 1996)
# GABAa
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
Alpha_AMPA = 3.94
Beta_AMPA = 0.18
Erev_AMPA = 0.
Deadtime_AMPA = 1.
# ------------------------------------------------------------------------------

# HPC AXONAL DELAYS---------------------------------------------------
dx = 10e-3         # 10 microns = 10*10^-3 mm
V_ax_P = 0.5       # PY axonal conductance velocity (mm/msec)
V_ax_I = 0.1       # IN axonal conductance velocity
#-------------------------------------------------------------------------------

# CONNECTIVITY Parameters-------------------------------------------------------
kPPh = 200. 
kPIh = 20.   
kIPh = 400. 
kIIh = 100.
sigma_PYh = 100.   # PY connectivity deviation in cell numbers
sigma_INh = 30.    # IN connectivity deviation in cell numbers
#-------------------------------------------------------------------------------

# CA3 CONNECTIVITY Parameters---------------------------------------------------
# kPP3 = 55. 
# kPI3 = 5.   
# kIP3 = 68. 
# sigma_PY3 = 100.   # PY connectivity deviation in cell numbers
# sigma_IN3 = 30.    # IN connectivity deviation in cell numbers
# kPP3 = 200. 
# kPI3 = 20.   
# kIP3 = 400. 
# sigma_PY3 = 100.   # PY connectivity deviation in cell numbers
# sigma_IN3 = 30.    # IN connectivity deviation in cell numbers
#-------------------------------------------------------------------------------

# # CA1 CONNECTIVITY Parameters---------------------------------------------------
# kPI1 = 20.  
# kIP1 = 400.
# kII1 = 100.
# sigma_PY1 = 100.   # PY connectivity deviation in cell numbers
# sigma_IN1 = 10.    # IN connectivity deviation in cell numbers
# #-------------------------------------------------------------------------------

# # SCHAFFER CONNECTIVITY Parameters----------------------------------------------
# k31 = 130.         # number of connections per CA3 PY cell
# sigma_31 = 120.    # connectivity deviation in cell numbers
# multiP = 20.       # number of multiple synapses on each connection
# multiI = 20.
#-------------------------------------------------------------------------------