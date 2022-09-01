from brian2 import *

# mV, cm, ms, msiemens/cm**2, uA/cm**2, ufarad/cm**2, mmole   

# NETWORK & CELL SIZES (in cm**2)-----------------------------------------------
N_PYc = 1000        
N_INc = 250


Gc = N_PYc + N_INc     # total number of cells
L = 0.5                # total length of PFC network (5*mm)

As = 1.5e-4            # Dendrtic compartment size (0.015*(mm**2))
Ad = 3.5e-4            # Axosomatic compartment size (0.035*(mm**2)) 
S_PY = As + Ad         # Total pyramidal cell size (0.05*(mm**2))
S_IN = 2e-4            # Total interneuron size (0.02*(mm**2))
#-------------------------------------------------------------------------------

# CORTICAL PYRAMIDAL CELLS Parameters-------------------------------------------
C_c = 1.   
phi = 4.      
gNa_c = 50. 
gK_c = 10.5  
gA_c = 1.   
gKS_c = 0.576  
gKNa_c = 1.33     #1.33 #0.5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! con 1.33 meno spike
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
# gAMPA_PP_c = 5.4*(1.e-6/Ad)    
# gAMPA_PI_c = 2.25*(1.e-6/S_IN)
# gGABA_IP_c = 4.15*(1.e-6/As)  
# gGABA_II_c = 0.165*(1.e-6/S_IN)

# gAMPA_PP_c = 50* 5.4*(1.e-6/Ad)    
# gAMPA_PI_c = 20*2.25*(1.e-6/S_IN)
# gGABA_IP_c = 4.15*(1.e-6/As)  
# gGABA_II_c = 20*0.165*(1.e-6/S_IN)

# gAMPA_PP_c = 30* 5.4*(1.e-6/Ad)    #29.3, 29.5, 30
# gAMPA_PI_c = 20*2.25*(1.e-6/S_IN)
# gGABA_IP_c = 4.15*(1.e-6/As)  #1-6
# gGABA_II_c = 5*0.165*(1.e-6/S_IN) 

gAMPA_PP_c = 30* 5.4*(1.e-6/Ad)    #29.3, 29.5, 30
gAMPA_PI_c = 10*2.25*(1.e-6/S_IN)
gGABA_IP_c = 4.15*(1.e-6/As)  #1-6
gGABA_II_c = 10*0.165*(1.e-6/S_IN) 

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
Alpha_AMPA = 1.54
# Alpha_AMPA = 5.54
Beta_AMPA = 0.18
Erev_AMPA = 0.
Deadtime_AMPA = 1.
# ------------------------------------------------------------------------------

# PFC CONNECTIVITY Parameters---------------------------------------------------
kPPc = kPIc = kIPc = kIIc = 20. 
sigma_PYc = 0.0250/(L/Gc)  # PY connectivity deviation in cell numbers
sigma_INc = 0.0125/(L/Gc)  # IN connectivity deviation in cell numbers
#-------------------------------------------------------------------------------