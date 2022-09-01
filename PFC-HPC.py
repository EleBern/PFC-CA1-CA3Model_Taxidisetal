from brian2 import *
from brian2tools import *
from numpy import random
import scipy.io
from scipy.io import savemat
import time
seed(6)

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
# refractoriness for HH type cells. Refractory was 3*ms for PYh and 1*ms for IN.
PYc = NeuronGroup(N_PYc, eqs_PYc, threshold= 'Vs>0', refractory='Vs>=0', order=2, method = 'euler')
PYh = NeuronGroup(N_PYh, eqs_PYh, threshold='Vs>0', refractory='Vs>=0', order=2, method = 'euler')
IN = NeuronGroup(N_IN, eqs_IN, threshold='Vi>0', refractory='Vi>0', order=2, method = 'euler')

INc = IN[0:N_INc]        # PFC INs
INh = IN[N_INc:N_IN]     # HPC INs
#-------------------------------------------------------------------------------

# VARYING INTRINSIC CELL PARAMETERS --------------------------------------------
PYc.VL = VL_c + 0.3*randn(N_PYc)      
PYc.gL = gL_c + 0.0067*randn(N_PYc)  
PYc.gsd = gsd_c + 1.e-4*randn(N_PYc)     #*msiemens

PYh.VL = VL_h + (VL_h/200.)*randn(N_PYh)
PYh.gL = gL_h + (gL_h/200.)*randn(N_PYh)
PYh.gsd = gsd_h + (gsd_h/200.)*randn(N_PYh)

IN.VLi = VLi_m + (VLi_m/200.)*randn(N_IN)     
IN.gLi = gLi_m + (gLi_m/200.)*randn(N_IN)  

PYh.gCa_h = 'int(i%2==1)*gCa_1 + int(i%2==0)*gCa_3' # odd neurons CA1, even neurons CA3
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

# PY3 Initialization - potrei farla unica x entrambe le aree
PYh.Vs = 'int(i%2==0)+(-65.60479241 + 6.6*randn())'
PYh.Vd = 'int(i%2==0)+(-65.73747962 + 6.6*randn())'
PYh.Ca = 'int(i%2==0)+(0.84537446 + 0.08*randn())'
PYh.h = 'int(i%2==0)+(0.99903134 + 0.1*randn())'
PYh.n = 'int(i%2==0)+(0.00039902 + 0.00004*randn())'
PYh.s = 'int(i%2==0)+(0.00851051 + 0.0008*randn())'
PYh.c = 'int(i%2==0)+(0.00630044 + 0.0006*randn())'
PYh.q = 'int(i%2==0)+(0.163881 + 0.02*randn())'

# PY1 Initialization
PYh.Vs = 'int(i%2==1)+(-62.89223689 + 6.3*randn())'
PYh.Vd = 'int(i%2==1)+(-62.98248752 + 6.3*randn())'
PYh.Ca = 'int(i%2==1)+(0.21664282 + 0.02*randn())'
PYh.h = 'int(i%2==1)+(0.99806345 + 0.1*randn())'
PYh.n = 'int(i%2==1)+(0.00068604 + 0.00007*randn())'
PYh.s = 'int(i%2==1)+(0.01086703 + 0.001*randn())'
PYh.c = 'int(i%2==1)+(0.00809387 + 0.0008*randn())'
PYh.q = 'int(i%2==1)+(0.0811213 + 0.008*randn())'

IN.Vi = -61.95897754 + 6.2*randn(N_IN)
IN.hi = 0.72463045 + 0.07*randn(N_IN)
IN.ni = 0.10410982 + 0.01*randn(N_IN)
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
CPPc.Alpha_AMPA = Alpha_AMPA_c

sourcesCPIc, targetsCPIc = CPI['CPI'].nonzero()
CPIc = Synapses(PYc, INc, model = SynAMPA_PI, on_pre=AMPA_PI_fire, method = 'euler')
CPIc.connect(i=sourcesCPIc, j=targetsCPIc)
CPIc.gAMPA_PI = 'gAMPA_PI_c / N_incoming'
CPIc.Alpha_AMPA = Alpha_AMPA_c

sourcesCIPc, targetsCIPc = CIP['CIP'].nonzero()
CIPc = Synapses(INc, PYc, model = SynGABA_IP, on_pre= GABA_IP_fire, method = 'euler')
CIPc.connect(i=sourcesCIPc, j=targetsCIPc)
CIPc.gGABA_IP = 'gGABA_IP_c / N_incoming'

sourcesCIIc, targetsCIIc = CII['CII'].nonzero()
CIIc = Synapses(INc, INc, model = SynGABA_II, on_pre= GABA_II_fire, method = 'euler')
CIIc.connect(i=sourcesCIIc, j=targetsCIIc)
CIIc.gGABA_II = 'gGABA_II_c / N_incoming'
#-------------------------------------------------------------------------------

# HPC CONNECTIVITY--------------------------------------------------------------
# from Connectivity_HP import C_PP_h,C_PI_h,C_II_h,C_IP_h,d_PP_h,d_II_h,d_PI_h     # import the connectivity matrices
C_PP_h = scipy.io.loadmat('C_PP_h.mat')
C_PI_h = scipy.io.loadmat('C_PI_h.mat')
C_II_h = scipy.io.loadmat('C_II_h.mat')
C_IP_h = scipy.io.loadmat('C_IP_h.mat')

d_PP_h = scipy.io.loadmat('d_PP_h.mat')
d_II_h = scipy.io.loadmat('d_II_h.mat')
d_PI_h = scipy.io.loadmat('d_PI_h.mat')
d_IP_h = scipy.io.loadmat('d_IP_h.mat')

sourcesCPPh, targetsCPPh = C_PP_h['C_PP_h'].nonzero()
CPPh = Synapses(PYh, PYh, model = SynAMPA_PP, on_pre=AMPA_PP_fire, method = 'euler')
CPPh.connect(i=sourcesCPPh, j=targetsCPPh)
delayCPPh = abs(d_PP_h['d_PP_h'][sourcesCPPh,targetsCPPh])*dx/V_ax_P*ms
delayCPPh = np.where(delayCPPh < sigma_PYh*4*dx/V_ax_P*ms, delayCPPh, sigma_PYh*4*dx/V_ax_P*ms)
CPPh.delay = delayCPPh
         
sourcesCPIh, targetsCPIh = C_PI_h['C_PI_h'].nonzero()
CPIh = Synapses(PYh, INh, model = SynAMPA_PI, on_pre=AMPA_PI_fire, method = 'euler')
CPIh.connect(i=sourcesCPIh, j=targetsCPIh)
delayCPIh = abs(d_PI_h['d_PI_h'][sourcesCPIh,targetsCPIh])*dx/V_ax_P*ms
delayCPIh = np.where(delayCPIh < sigma_PYh*4*dx/V_ax_P*ms, delayCPIh, sigma_PYh*4*dx/V_ax_P*ms)
CPIh.delay = delayCPIh
               
sourcesCIPh, targetsCIPh = C_IP_h['C_IP_h'].nonzero()
CIPh = Synapses(INh, PYh, model = SynGABA_IP, on_pre=GABA_IP_fire, method = 'euler')
CIPh.connect(i=sourcesCIPh, j=targetsCIPh)
delayCIPh = abs(d_IP_h['d_IP_h'][sourcesCIPh,targetsCIPh])*dx/V_ax_I*ms

sourcesCIIh, targetsCIIh = C_II_h['C_II_h'].nonzero()
CIIh = Synapses(INh, INh, model = SynGABA_II, on_pre=GABA_II_fire, method = 'euler')
CIIh.connect(i=sourcesCIIh, j=targetsCIIh)
delayCIIh = abs(d_II_h['d_II_h'][sourcesCIIh,targetsCIIh])*dx/V_ax_I*ms
delayCIIh = np.where(delayCIIh < sigma_INh*5*dx/V_ax_I*ms, delayCIIh, sigma_INh*5*dx/V_ax_I*ms)
CIIh.delay = delayCIIh

CPPh.Alpha_AMPA = Alpha_AMPA_h
CPIh.Alpha_AMPA = Alpha_AMPA_h
CPPh.gAMPA_PP = 'gAMPA_PP_h / N_incoming'
CIPh.gGABA_IP = 'gGABA_IP_h / N_incoming' 
CPIh.gAMPA_PI = 'gAMPA_PI_h / N_incoming' 
CIIh.gGABA_II = 'gGABA_II_h / N_incoming'            
#-------------------------------------------------------------------------------

# PFC-HPC CONNECTIVITY----------------------------------------------------------
# from Connectivity_PFC_HPC import C_ch_P,C_ch_I
C_ch_P = scipy.io.loadmat('C_ch_P.mat')
C_ch_I = scipy.io.loadmat('C_ch_I.mat')

sourcesCch_P, targetsCch_P = C_ch_P['C_ch_P'].nonzero()
Cch_P = Synapses(PYc, PYh, model = SynC_HP_PP, on_pre=C_HP_PP_fire, method = 'euler')
Cch_P.connect(i=sourcesCch_P, j=targetsCch_P)
Cch_P.gAMPA_PP = 'gAMPA_PP_h / N_incoming'
Cch_P.Alpha_AMPA = Alpha_AMPA_h

sourcesCch_I, targetsCch_I = C_ch_I['C_ch_I'].nonzero()
Cch_I = Synapses(PYc, INh, model = SynC_HP_PI, on_pre=C_HP_PI_fire, method = 'euler')
Cch_I.connect(i=sourcesCch_I, j=targetsCch_I)
Cch_I.gAMPA_PI = 'gAMPA_PI_h / N_incoming' 
Cch_I.Alpha_AMPA = Alpha_AMPA_h
#-------------------------------------------------------------------------------

# HPC-PFC CONNECTIVITY----------------------------------------------------------
# from Connectivity_HPC_PFC import C_hc_P,C_hc_I
C_hc_P = scipy.io.loadmat('C_hc_P.mat')
C_hc_I = scipy.io.loadmat('C_hc_I.mat')

sourcesChc_P, targetsChc_P = C_hc_P['C_hc_P'].nonzero()
Chc_P = Synapses(PYh, PYc, model = SynC_HP_PP, on_pre=C_HP_PP_fire, method = 'euler')
Chc_P.connect(i=sourcesChc_P, j=targetsChc_P)
Chc_P.gAMPA_PP = 'gAMPA_PP_c / N_incoming'
Chc_P.Alpha_AMPA = Alpha_AMPA_c

sourcesChc_I, targetsChc_I = C_hc_I['C_hc_I'].nonzero()
Chc_I = Synapses(PYh, INc, model = SynC_HP_PI, on_pre=C_HP_PI_fire, method = 'euler')
Chc_I.connect(i=sourcesChc_I, j=targetsChc_I)
Chc_I.gAMPA_PI = 'gAMPA_PI_c / N_incoming' 
Chc_I.Alpha_AMPA = Alpha_AMPA_c
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
spikes_ph = SpikeMonitor(PYh,record=True)
spikes_ih = SpikeMonitor(INh,record=True)

V_pc = StateMonitor(PYc,'Vs',record=[0, 200, 400, 600, 800, 999], dt=0.05*ms)
V_ic = StateMonitor(INc,'Vi',record=[0, 50, 100, 150, 200, 249], dt=0.05*ms)
V_p3 = StateMonitor(PYh,'Vs',record=[0, 400, 800, 1200, 1600, 1998], dt=0.05*ms)
V_ih = StateMonitor(IN,'Vi',record=[0, 50, 100, 150, 199], dt=0.05*ms)
V_p1 = StateMonitor(PYh,'Vs',record=[1, 401, 801, 1201, 1601, 1999], dt=0.05*ms)
#-------------------------------------------------------------------------------

# NOISE-------------------------------------------------------------------------
PYc.run_regularly('Id = randn()*0.7/(sim_clock/ms)**0.5')
PYc.run_regularly('Is = randn()*0.7/(sim_clock/ms)**0.5')
INc.run_regularly('Ii = randn()*0.3/(sim_clock/ms)**0.5')

PYh.run_regularly('Id = randn()*0.2/(sim_clock/ms)**0.5')
PYh.run_regularly('Is = randn()*0.2/(sim_clock/ms)**0.5')
INh.run_regularly('Ii = randn()*0.2/(sim_clock/ms)**0.5')
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
CPPh.active = False
CPIh.active = False
CIPh.active = False
CIIh.active = False
Cch_P.active = False
Cch_I.active = False
Chc_P.active = False
Chc_I.active = False


spikes_pc.active = False
spikes_ic.active = False
V_pc.active = False
V_ic.active = False
spikes_ph.active = False
spikes_ih.active = False
V_p3.active = False
V_p1.active = False
V_ih.active = False

print_time.active = False

net.run(1*second)

# Connect neurons
print ("Connecting neurons...")
CPPc.active = True
CPIc.active = True
CIPc.active = True
CIIc.active = True
CPPh.active = True
CPIh.active = True
CIPh.active = True
CIIh.active = True
Cch_P.active = True
Cch_I.active = True
Chc_P.active = True
Chc_I.active = True

net.run(1*second)

# Monitor V & spikes after V has stabilized after making connections
spikes_pc.active = True
spikes_ic.active = True
V_pc.active = True
V_ic.active = True
spikes_ph.active = True
spikes_ih.active = True
V_p3.active = True
V_p1.active = True
V_ih.active = True

print_time.active = True

start_time=time.time()
net.run(t_final)
duration=time.time()-start_time
print ("Simulation time:",duration,"seconds")
#-------------------------------------------------------------------------------

# EXPORT DATA-------------------------------------------------------------------
savemat('V_pc.mat', {'t':array(V_pc.t), 'Vs':array(V_pc.Vs)})
savemat('V_ic.mat', {'t':array(V_ic.t), 'Vi':array(V_ic.Vi)})
savemat('V_p3.mat', {'t':array(V_p3.t), 'Vs':array(V_p3.Vs)})
savemat('V_p1.mat', {'t':array(V_p1.t), 'Vs':array(V_p1.Vs)})
savemat('V_ih.mat', {'t':array(V_ih.t), 'Vi':array(V_ih.Vi)})

savemat('Spikes_pc.mat', {'Neuron_index':array(spikes_pc.i), 'Spike_time':array( spikes_pc.t)})
savemat('Spikes_ic.mat', {'Neuron_index':array(spikes_ic.i), 'Spike_time':array(spikes_ic.t)})
savemat('Spikes_ih.mat', {'Neuron_index':array(spikes_ih.i), 'Spike_time':array(spikes_ih.t)})
savemat('Spikes_ph.mat', {'Neuron_index':array(spikes_ph.i), 'Spike_time':array(spikes_ph.t)})
#-------------------------------------------------------------------------------

# PLOTS-------------------------------------------------------------------------
# Raster plot pyramidal cells cortex
figure()
title("Spikes PY cortex")
brian_plot(spikes_pc)

# (Vs,t) pyramidal cells cortex
figure()
title("V PY cortex")
brian_plot(V_pc[0])
figure()
title("V PY cortex")
brian_plot(V_pc[600])
figure()
title("V PY cortex")
brian_plot(V_pc[999])

# Raster plot interneurons cortex
figure()
title("Spikes IN cortex")
brian_plot(spikes_ic)

# (Vs,t) interneurons cortex
figure()
title("V IN cortex")
brian_plot(V_ic[0])
figure()
title("V IN cortex")
brian_plot(V_ic[100])
figure()
title("V IN cortex")
brian_plot(V_ic[200])

# Raster plot pyramidal cells HPC
figure()
title("Spikes PY HPC")
brian_plot(spikes_ph)

# (Vs,t) pyramidal cells CA1
figure()
title("V PY CA1")
brian_plot(V_p1[1])
figure()
title("V PY CA1")
brian_plot(V_p1[801])
figure()
title("V PY CA1")
brian_plot(V_p1[1999])

# (Vs,t) pyramidal cells CA3
figure()
title("V PY CA3")
brian_plot(V_p3[0])
figure()
title("V PY CA3")
brian_plot(V_p3[800])
figure()
title("V PY CA3")
brian_plot(V_p3[1998])

# Raster plot interneurons HPC
figure()
title("Spikes IN HPC")
brian_plot(spikes_ih)

# (Vs,t) interneurons HPC
figure()
title("V IN CA3")
brian_plot(V_ih[0])
figure()
title("V IN CA3")
brian_plot(V_ih[100])
figure()
title("V IN CA3")
brian_plot(V_ih[199])
#-------------------------------------------------------------------------------

# Stats-------------------------------------------------------------------------
# Cortex
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

MFR_totalc = (MFR_pc[N_PYc] + MFR_ic[N_INc]) / ((N_PYc + N_INc) * t_final)
MFR_total = MFR_pc[N_PYc] + MFR_ic[N_INc]
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

# Hippocampus
# Pyramidal cells CA1
spike_dict_p1 = spikes_ph.all_values()
[spike_dict_p1['t'].pop(key) for key in np.arange(1,2000,2)]
NofSpikes_p1 = [0]*(N_PY1)
MFR_p1 = [0]*(N_PY1+1)
Var_p1 = 0.
active_p1 = 0

# Pyramidal cells CA3
spike_dict_p3 = spikes_ph.all_values()
[spike_dict_p3['t'].pop(key) for key in np.arange(0,2000,2)]
NofSpikes_p3 = [0]*(N_PY3)
MFR_p3 = [0]*(N_PY3+1)
Var_p3 = 0.
active_p3 = 0

# Interneurons
spike_dict_ih = spikes_ih.all_values()
NofSpikes_ih = [0]*(N_INh)
MFR_ih = [0]*(N_INh+1)
Var_ih = 0.
active_ih = 0

# Mean firing rate
# Pyramidal cells CA1
x = 0
for x in range(N_PY1):
    NofSpikes_p1[x] = len(spike_dict_p1['t'][x*2])
    MFR_p1[x] = NofSpikes_p1[x]/t_final
    MFR_p1[N_PY1] = MFR_p1[N_PY1] + NofSpikes_p1[x]
    if MFR_p1[x] > 0:
        active_p1 = active_p1 + 1
        
# Pyramidal cells CA3
x = 0
for x in range(N_PY3):
    NofSpikes_p3[x] = len(spike_dict_p3['t'][x*2+1])
    MFR_p3[x] = NofSpikes_p3[x]/t_final
    MFR_p3[N_PY3] = MFR_p3[N_PY3] + NofSpikes_p3[x]
    if MFR_p3[x] > 0:
        active_p3 = active_p3 + 1

# Interneurons
x = 0
for x in range(N_INh):
    NofSpikes_ih[x] = len(spike_dict_ih['t'][x])
    MFR_ih[x] = NofSpikes_ih[x]/t_final
    MFR_ih[N_INh] = MFR_ih[N_INh] + NofSpikes_ih[x]
    if MFR_ih[x] > 0:
        active_ih = active_ih + 1

MFR_totalh = (MFR_p1[N_PY1] + MFR_p3[N_PY3] + MFR_ih[N_INh]) / ((N_PYh + N_INh) * t_final)
MFR_total = (MFR_total + MFR_p1[N_PY1] + MFR_p3[N_PY3] + MFR_ih[N_INh]) / ((N_PYc + N_PYh + N_IN) * t_final)
MFR_p1[N_PY1] = MFR_p1[N_PY1]/(N_PY1*t_final)  
MFR_p3[N_PY3] = MFR_p3[N_PY3]/(N_PY3*t_final)  
MFR_ih[N_INh] = MFR_ih[N_INh]/(N_INh*t_final) 

# Variance 
# Pyramidal cells CA1
x = 0
for x in range(N_PY1):
    Var_p1 = Var_p1 + (MFR_p1[x] - MFR_p1[N_PY1])**2

# Pyramidal cells CA3
x = 0
for x in range(N_PY3):
    Var_p3 = Var_p3 + (MFR_p3[x] - MFR_p3[N_PY3])**2

# Interneurons
x = 0
for x in range(N_INh):
    Var_ih = Var_ih + (MFR_ih[x] - MFR_ih[N_INh])**2
    
Var_p1 = Var_p1 / N_PY1
Var_p3 = Var_p3 / N_PY3
Var_ih = Var_ih / N_INh

print('MFR of PYc cells: ', MFR_pc[N_PYc])
print('Variance of PYc cells: ', Var_pc)
print('MFR of only active PYc cells: ', MFR_pc[N_PYc]*N_PYc/active_pc)
print('Percentage of active PYc cells: ', active_pc*100/N_PYc, ' %')
print('MFR of INc cells: ', MFR_ic[N_INc])
print('Variance of INc cells: ', Var_ic)
print('MFR of only active INc cells: ', MFR_ic[N_INc]*N_INc/active_ic)
print('Percentage of active INc cells: ', active_ic*100/N_INc,' %')
print('MFR of all cortex cells: ', MFR_totalc)
print('Percentage of total cortical active cells: ', (active_ic+active_pc)*100/(N_INc+N_PYc),' %')  

print('MFR of PY1 cells: ', MFR_p1[N_PY1])
print('Variance of PY1 cells: ', Var_p1)
print('MFR of PY3 cells: ', MFR_p3[N_PY3])
print('Variance of PY3 cells: ', Var_p3)
print('MFR of INh cells: ', MFR_ih[N_INh])
print('Variance of INh cells: ', Var_ih)
print('MFR of all hippocampal cells: ', MFR_totalh)
print('MFR of only active PY1 cells: ', MFR_p1[N_PY1]*N_PY1/active_p1)
print('MFR of only active PY3 cells: ', MFR_p3[N_PY3]*N_PY3/active_p3)
print('MFR of only active INh cells: ', MFR_ih[N_INh]*N_INh/active_ih)
print('Percentage of active PY1 cells: ', active_p1*100/N_PY1, ' %')
print('Percentage of active PY3 cells: ', active_p3*100/N_PY3, ' %')
print('Percentage of active PY cells: ', (active_p1+active_p3)*100/(N_PYh),' %')
print('Percentage of active IN cells: ', active_ih*100/(N_INh),' %')
print('Percentage of total hippocampal active cells: ', (active_ih+active_p1+active_p3)*100/(N_INh+N_PYh),' %')  

print('MFR of all cells: ', MFR_total)   
print('Percentage of total active cells: ', (active_ic+active_pc+active_ih+active_p1+active_p3)*100/(N_IN+N_PYc+N_PYh),' %') 
#-------------------------------------------------------------------------------  
