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
# refractoriness for HH type cells. Refractory was 3*ms for PYh and 1*ms for IN.
PYh = NeuronGroup(N_PYh, eqs_PYh, threshold='Vs > 0', refractory='Vs >= 0', order=2, method = 'euler')
INh = NeuronGroup(N_IN, eqs_IN, threshold='Vi > 0', refractory='Vi >= 0', order=2, method = 'euler')
#-------------------------------------------------------------------------------

# VARYING INTRINSIC CELL PARAMETERS --------------------------------------------
PYh.VL = VL_h + (VL_h/200.)*randn(N_PYh)
PYh.gL = gL_h + (gL_h/200.)*randn(N_PYh)
PYh.gsd = gsd_h + (gsd_h/200.)*randn(N_PYh)

INh.VLi = VLi_m + (VLi_m/200.)*randn(N_IN)     
INh.gLi = gLi_m + (gLi_m/200.)*randn(N_IN)  

PYh.gCa_h = 'int(i%2==1)*gCa_1 + int(i%2==0)*gCa_3' # odd neurons CA1, even neurons CA3
#-------------------------------------------------------------------------------

# INITIALIZATION----------------------------------------------------------------
# From Taidis

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

# Interneurons
INh.Vi = -61.95897754 + 6.2*randn(N_IN)
INh.hi = 0.72463045 + 0.07*randn(N_IN)
INh.ni = 0.10410982 + 0.01*randn(N_IN)
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

CPPh.gAMPA_PP = 'gAMPA_PP_h / N_incoming'
CIPh.gGABA_IP = 'gGABA_IP_h / N_incoming' 
CPIh.gAMPA_PI = 'gAMPA_PI_h / N_incoming' 
CIIh.gGABA_II = 'gGABA_II_h / N_incoming'                  
#-------------------------------------------------------------------------------

# NETWORK OPERATIONS--------------------------------------------------------
@network_operation(dt=print_clock)     # print every second of simulation time
def print_time():
    global to_print
    to_print += print_clock
    print ("Simulating second ", to_print)
#-------------------------------------------------------------------------------

# MONITORS----------------------------------------------------------------------
spikes_p = SpikeMonitor(PYh,record=True)
spikes_i = SpikeMonitor(INh,record=True)

V_p3 = StateMonitor(PYh,'Vs',record=[0, 400, 800, 1200, 1600, 1998], dt=0.05*ms)
V_i = StateMonitor(INh,'Vi',record=[0, 50, 100, 150, 199], dt=0.05*ms)
V_p1 = StateMonitor(PYh,'Vs',record=[1, 401, 801, 1201, 1601, 1999], dt=0.05*ms)
#-------------------------------------------------------------------------------

# NOISE-------------------------------------------------------------------------
PYh.run_regularly('Id = randn()*0.2/(sim_clock/ms)**0.5')
PYh.run_regularly('Is = randn()*0.2/(sim_clock/ms)**0.5')
INh.run_regularly('Ii = randn()*0.2/(sim_clock/ms)**0.5')
#-------------------------------------------------------------------------------

# RUN---------------------------------------------------------------------------
print ("Network construction time:",time.time()-start_time,"seconds")
print ("Simulation running...")

net = Network(collect())

#  Wait for V to stabilize before connecting neurons 
CPPh.active = False
CPIh.active = False
CIPh.active = False
CIIh.active = False
spikes_p.active = False
spikes_i.active = False
V_p3.active = False
V_p1.active = False
V_i.active = False
print_time.active = False

net.run(1*second)

# Connect neurons
print ("Connecting neurons...")
CPPh.active = True
CPIh.active = True
CIPh.active = True
CIIh.active = True

net.run(1*second)

# Monitor V & spikes after V has stabilized after making connections
spikes_p.active = True
spikes_i.active = True
V_p3.active = True
V_p1.active = True
V_i.active = True

print_time.active = True

start_time=time.time()
net.run(t_final)
duration=time.time()-start_time
print ("Simulation time:",duration,"seconds")
#-------------------------------------------------------------------------------

# EXPORT DATA-------------------------------------------------------------------
savemat('V_p3.mat', {'t':array(V_p3.t), 'Vs':array(V_p3.Vs)})
savemat('V_p1.mat', {'t':array(V_p1.t), 'Vs':array(V_p1.Vs)})
savemat('V_i.mat', {'t':array(V_i.t), 'Vi':array(V_i.Vi)})

savemat('Spikes_i.mat', {'Neuron_index':array(spikes_i.i), 'Spike_time':array(spikes_i.t)})
savemat('Spikes_p.mat', {'Neuron_index':array(spikes_p.i), 'Spike_time':array(spikes_p.t)})
#-------------------------------------------------------------------------------

# PLOTS-------------------------------------------------------------------------
# Raster plot pyramidal cells
figure()
title("Spikes PY HPC")
brian_plot(spikes_p)

# Raster plot interneurons
figure()
title("Spikes IN HPC")
brian_plot(spikes_i)

# (Vs,t) pyramidal cells CA1
figure()
title("V PY CA1")
brian_plot(V_p1[1])
figure()
title("V PY CA1")
brian_plot(V_p1[401])
figure()
title("V PY CA1")
brian_plot(V_p1[1999])

# (Vs,t) interneurons HPC
figure()
title("V IN")
brian_plot(V_i[0])
figure()
title("V IN")
brian_plot(V_i[50])
figure()
title("V IN")
brian_plot(V_i[199])

# (Vs,t) pyramidal cells CA3
figure()
title("V PY CA3")
brian_plot(V_p3[0])
figure()
title("V PY CA3")
brian_plot(V_p3[400])
figure()
title("V PY CA3")
brian_plot(V_p3[1998])

# Stats-------------------------------------------------------------------------
# Pyramidal cells CA1
spike_dict_p1 = spikes_p.all_values()
[spike_dict_p1['t'].pop(key) for key in np.arange(1,2000,2)]
NofSpikes_p1 = [0]*(N_PY1)
MFR_p1 = [0]*(N_PY1+1)
Var_p1 = 0.
active_p1 = 0

# Pyramidal cells CA3
spike_dict_p3 = spikes_p.all_values()
[spike_dict_p3['t'].pop(key) for key in np.arange(2,2000,2)]
NofSpikes_p3 = [0]*(N_PY3)
MFR_p3 = [0]*(N_PY3+1)
Var_p3 = 0.
active_p3 = 0

# Interneurons
spike_dict_ih = spikes_i.all_values()
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

MFR_total = (MFR_p1[N_PY1] + MFR_p3[N_PY3] + MFR_ih[N_INh]) / ((N_PYh + N_INh) * t_final)
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

print('MFR of PY1 cells: ', MFR_p1[N_PY1])
print('Variance of PY1 cells: ', Var_p1)
print('MFR of PY3 cells: ', MFR_p3[N_PY3])
print('Variance of PY3 cells: ', Var_p3)
print('MFR of INh cells: ', MFR_ih[N_INh])
print('Variance of INh cells: ', Var_ih)
print('MFR of all cells: ', MFR_total)
print('MFR of only active PY1 cells: ', MFR_p1[N_PY1]*N_PY1/active_p1)
print('MFR of only active PY3 cells: ', MFR_p3[N_PY3]*N_PY3/active_p3)
print('MFR of only active INh cells: ', MFR_ih[N_INh]*N_INh/active_ih)
print('Percentage of active PY1 cells: ', active_p1*100/N_PY1, ' %')
print('Percentage of active PY3 cells: ', active_p3*100/N_PY3, ' %')
print('Percentage of active PY cells: ', (active_p1+active_p3)*100/(N_PYh),' %')
print('Percentage of active IN cells: ', active_ih*100/(N_INh),' %')
print('Percentage of total active cells: ', (active_ih+active_p1+active_p3)*100/(N_INh+N_PYh),' %')  
#-------------------------------------------------------------------------------  