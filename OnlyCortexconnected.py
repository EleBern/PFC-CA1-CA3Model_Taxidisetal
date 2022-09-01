# -*- coding: utf-8 -*-
"""
Created on Tue May  4 14:08:13 2021

@author: ely97
"""

from brian2 import *
from brian2tools import *
from numpy import random
import scipy.io
from scipy.io import savemat
import time
seed(20)

from Parameters import *
from Equations import *       

prefs.codegen.target = "numpy"

# TIME Parameters---------------------------------------------------------------
sim_clock = 0.01*ms         # simulation timestep 
defaultclock.dt = sim_clock # default simulation timestep for all objects!!
rec_clock = 1*ms            # monitors recording timestep
print_clock = 1*second      # printing-seconds clock
to_print = 0*second
start_time = time.time()    # start recording real time
t_final = 60*second          # simulation duration
#-------------------------------------------------------------------------------

# GROUPS------------------------------------------------------------------------
# Changed threshold so it's as in Destexhe model and refractory as in Brian2 user guide
# refractoriness for HH type cells.
PYc = NeuronGroup(N_PYc, eqs_PYc, threshold= 'Vs > 0', refractory= 'Vs >= 0', order=2, method='euler')
# In PYc Taxidis - threshold = 'Vs > 10'. Cambiato a 0 seguendo parametri Destexhe.
INc = NeuronGroup(N_INc, eqs_IN, threshold='Vi > 0', refractory= 'Vi >= 0', order=2, method='euler')
# In INc Taxidis - threshold = 'Vi > -10'. Cambiato a 0 seguendo parametri Destexhe.
#-------------------------------------------------------------------------------

# VARYING INTRINSIC CELL PARAMETERS --------------------------------------------
PYc.VL = VL_c + 0.3*randn(N_PYc)      
PYc.gL = gL_c + 0.0067*randn(N_PYc)  
PYc.gsd = gsd_c + 1.e-4*randn(N_PYc)     #*msiemens

# TAXIDIS
INc.VL = VLi_m + (VLi_m/200.)*randn(N_INc)     
INc.gL = gLi_m + (gLi_m/200.)*randn(N_INc)  
#-------------------------------------------------------------------------------

# INITIALIZATION----------------------------------------------------------------
PYc.Vd = -75.53641681 + 7.5*randn(N_PYc)
PYc.Vs = -75.47003318 + 7.5*randn(N_PYc)
PYc.Na = 9.51153476 + 0.9*randn(N_PYc)
PYc.Ca = 0. 
PYc.h = 0.9956734 + 0.1*randn(N_PYc)
PYc.n = 0.01491113 + 0.001*randn(N_PYc)
PYc.hA = 0.31973314 + 0.03*randn(N_PYc)
PYc.mKS = 0.00169225 + 0.0001*randn(N_PYc)


INc.Vi = -61.95897754 + 6.2*randn(N_INc)
INc.hi = 0.72463045 + 0.07*randn(N_INc)
INc.ni = 0.10410982 + 0.01*randn(N_INc) 
#-------------------------------------------------------------------------------

# PFC CONNECTIVITY--------------------------------------------------------------
# from Connectivity_PFC import CPP,CPI,CIP,CII      # import the connectivity matrices
CPP = scipy.io.loadmat('CPP.mat')
CPI = scipy.io.loadmat('CPI.mat')
CII = scipy.io.loadmat('CII.mat')
CIP = scipy.io.loadmat('CIP.mat')

sourcesCPPc, targetsCPPc = CPP['CPP'].nonzero()
CPPc = Synapses(PYc, PYc, model = SynAMPA_PP, on_pre= AMPA_PP_fire, method = 'euler')
CPPc.connect(i=sourcesCPPc, j=targetsCPPc)
CPPc.gAMPA_PP = 'gAMPA_PP_c / N_incoming'

sourcesCPIc, targetsCPIc = CPI['CPI'].nonzero()
CPIc = Synapses(PYc, INc, model = SynAMPA_PI, on_pre=AMPA_PI_fire, method = 'euler')
CPIc.connect(i=sourcesCPIc, j=targetsCPIc)
CPIc.gAMPA_PI = 'gAMPA_PI_c / N_incoming'

sourcesCIPc, targetsCIPc = CIP['CIP'].nonzero()
CIPc = Synapses(INc, PYc, model = SynGABA_IP, on_pre= GABA_IP_fire, method = 'euler')
CIPc.connect(i=sourcesCIPc, j=targetsCIPc)
CIPc.gGABA_IP = 'gGABA_IP_c / N_incoming'

sourcesCIIc, targetsCIIc = CII['CII'].nonzero()
CIIc = Synapses(INc, INc, model = SynGABA_II, on_pre= GABA_II_fire, method = 'euler')
CIIc.connect(i=sourcesCIIc, j=targetsCIIc)
CIIc.gGABA_II = 'gGABA_II_c / N_incoming'
#-------------------------------------------------------------------------------

# NETWORK OPERATIONS--------------------------------------------------------
@network_operation(dt=print_clock)     # print every second of simulation time
def print_time():
    global to_print
    to_print += print_clock
    print ("Simulating second ", to_print)
#-------------------------------------------------------------------------------

# MONITORS----------------------------------------------------------------------
spikes_pc = SpikeMonitor(PYc,record=True)
spikes_ic = SpikeMonitor(INc,record=True)
V_pc = StateMonitor(PYc,'Vs',record=[1, 200, 400, 600, 800, 999], dt = 0.05*ms)
V_ic = StateMonitor(INc,'Vi',record=[1, 50, 100, 150, 200, 249], dt = 0.05*ms)
#-------------------------------------------------------------------------------

# NOISE-------------------------------------------------------------------------
PYc.run_regularly('Id = randn()*0.7/(sim_clock/ms)**0.5')
PYc.run_regularly('Is = randn()*0.7/(sim_clock/ms)**0.5')
INc.run_regularly('Ii = randn()*0.3/(sim_clock/ms)**0.5')

# PYc.run_regularly('Id = randn()*1.0/(sim_clock/ms)**0.5')
# PYc.run_regularly('Is = randn()*0.5/(sim_clock/ms)**0.5')
# INc.run_regularly('Ii = randn()*0.5/(sim_clock/ms)**0.5')
#-------------------------------------------------------------------------------

# RUN---------------------------------------------------------------------------
print ("Network construction time:",time.time()-start_time,"seconds")
print ("Simulation running...")

net = Network(collect())

#  Wait for V to stabilize before connecting neurons 
CPPc.active = False
CPIc.active = False
CIPc.active = False
CIIc.active = False
spikes_pc.active = False
spikes_ic.active = False
V_pc.active = False
V_ic.active = False
print_time.active = False

net.run(1*second)

# Connect neurons
print ("Connecting neurons...")
CPPc.active = True
CPIc.active = True
CIPc.active = True
CIIc.active = True

net.run(1*second)

# Monitor V & spikes after V has stabilized after making connections
spikes_pc.active = True
spikes_ic.active = True
V_pc.active = True
V_ic.active = True
print_time.active = True

start_time=time.time()
net.run(t_final)
duration=time.time()-start_time
print ("Simulation time:",duration,"seconds")
#-------------------------------------------------------------------------------

# EXPORT DATA-------------------------------------------------------------------
savemat('V_pc.mat', {'t':array(V_pc.t), 'Vs':array(V_pc.Vs)})
savemat('V_ic.mat', {'t':array(V_ic.t), 'Vi':array(V_ic.Vi)})
savemat('Spikes_pc.mat', {'Neuron_index':array(spikes_pc.i), 'Spike_time':array( spikes_pc.t)})
savemat('Spikes_ic.mat', {'Neuron_index':array(spikes_ic.i), 'Spike_time':array(spikes_ic.t)})
#-------------------------------------------------------------------------------

# PLOTS-------------------------------------------------------------------------
# Raster plot pyramidal cells cortex
figure()
title("Spikes PY cortex")
brian_plot(spikes_pc)

# (Vs,t) pyramidal cells cortex
figure()
title("V PY cortex")
brian_plot(V_pc[400])
figure()
title("V PY cortex")
brian_plot(V_pc[200])
figure()
title("V PY cortex")
brian_plot(V_pc[600])

# Raster plot interneurons cortex
figure()
title("Spikes IN cortex")
brian_plot(spikes_ic)

# (Vs,t) interneurons cortex
figure()
title("V IN cortex")
brian_plot(V_ic[1])
figure()
title("V IN cortex")
brian_plot(V_ic[100])
figure()
title("V IN cortex")
brian_plot(V_ic[200])
#-------------------------------------------------------------------------------

# Stats-------------------------------------------------------------------------
print('Total number of spikes of PYc cells')
print(spikes_pc.num_spikes)
print('Total number of spikes of INc cells')
print(spikes_ic.num_spikes)
# Pyramidal cells cortex
spike_dict_pc = spikes_pc.all_values()
NofSpikes_pc = [0]*(N_PYc)
MFR_pc = [0]*(N_PYc+1)
Var_pc = 0.
active_pc = 0
# Interneurons cortex
spike_dict_ic = spikes_ic.all_values()
NofSpikes_ic = [0]*(N_INc)
MFR_ic = [0]*(N_INc+1)
Var_ic = 0.
active_ic = 0

# Mean firing rate
# Pyramidal cells cortex
x = 0
for x in range(N_PYc):
    NofSpikes_pc[x] = len(spike_dict_pc['t'][x])
    MFR_pc[x] = NofSpikes_pc[x]/t_final
    MFR_pc[N_PYc] = MFR_pc[N_PYc] + NofSpikes_pc[x]
    if MFR_pc[x] > 0:
        active_pc = active_pc + 1
# Interneurons cortex    
x = 0
for x in range(N_INc):
    NofSpikes_ic[x] = len(spike_dict_ic['t'][x])
    MFR_ic[x] = NofSpikes_ic[x]/t_final
    MFR_ic[N_INc] = MFR_ic[N_INc] + NofSpikes_ic[x]
    if MFR_ic[x] > 0:
        active_ic = active_ic + 1

MFR_total = (MFR_pc[N_PYc] + MFR_ic[N_INc]) / ((N_PYc + N_INc) * t_final)
MFR_pc[N_PYc] = MFR_pc[N_PYc]/(N_PYc*t_final)  
MFR_ic[N_INc] = MFR_ic[N_INc]/(N_INc*t_final) 

# Variance 
# Pyramidal cells cortex
x = 0
for x in range(N_PYc):
    Var_pc = Var_pc + (MFR_pc[x] - MFR_pc[N_PYc])**2
# Interneurons cortex
x = 0
for x in range(N_INc):
    Var_ic = Var_ic + (MFR_ic[x] - MFR_ic[N_INc])**2
    
Var_pc = Var_pc / N_PYc
Var_ic = Var_ic / N_INc

print('MFR of PYc cells: ', MFR_pc[N_PYc])
print('Variance of PYc cells: ', Var_pc)
print('MFR of INc cells: ', MFR_ic[N_INc])
print('Variance of INc cells: ', Var_ic)
print('MFR of all cells: ', MFR_total)
print('MFR of only active PYc cells: ', MFR_pc[N_PYc]*N_PYc/active_pc)
print('MFR of only active INc cells: ', MFR_ic[N_INc]*N_INc/active_ic)
print('Percentage of active PYc cells: ', active_pc*100/N_PYc, ' %')
print('Percentage of active INc cells: ', active_ic*100/N_INc,' %')
print('Percentage of total active cells: ', (active_ic+active_pc)*100/(N_INc+N_PYc),' %') 
#-------------------------------------------------------------------------------