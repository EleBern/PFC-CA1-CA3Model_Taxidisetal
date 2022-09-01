from brian2 import *
from Parameters import *

# HPC PYRAMIDAL CELL Equations--------------------------------------------------
eqs_PYh = '''
    dVs/dt = (-( I_Ls + I_Na + I_K ) - I_sd/p - I_GABA_IP/p + Is/p)/C_h /ms      : 1
    dVd/dt = (-( I_Ld + I_Ca + I_KAHP + I_KCa ) + I_sd/(1.-p) - I_AMPA_PP/(1.-p) + Id/(1.-p))/C_h /ms : 1
    dCa/dt = (-0.13*I_Ca - 0.075*Ca) /ms                                         : 1
    
#SOMATIC currents    
# Leakage Current    
    I_Ls = gL*(Vs-VL)                                                            : 1

# Sodium Current  
    m_inf = alpham / (alpham + betam)                                            : 1
    dh/dt = (alphah - h*(alphah + betah)) /ms                                    : 1
    alpham = 0.32*4. / exprel((-46.9-Vs)/4.)                                     : 1
    betam = 0.28*5. / (exprel((Vs+19.9)/5.))                                     : 1 
    alphah = 0.128*exp((-43-Vs)/18.)                                             : 1
    betah = 4./(1+exp((-20-Vs)/5.))                                              : 1   
    I_Na = gNa_h*(m_inf**2)*h*(Vs-VNa_h)                                         : 1

# Potassium Current   
    dn/dt = (alphan - n*(alphan + betan)) /ms                                    : 1
    alphan = 0.016*5. / (exprel((-24.9-Vs)/5.))                                  : 1
    betan = 0.25*exp(-1-0.025*Vs)                                                : 1
    I_K = gK_h*n*(Vs-VK_h)                                                       : 1
    
  
#DENDRITIC Currents
# Leakage Current
    I_Ld = gL*(Vd-VL)                                                            : 1
    
# Ca+ Current
    ds/dt = (alphas - s*(alphas + betas)) /ms                                    : 1
    alphas = 1.6 / (exp(-0.072*(Vd-5))+1.)                                       : 1
    betas = 0.02*5. / (exprel((Vd+8.9)/5.))                                      : 1
    I_Ca = gCa_h*(s**2)*(Vd-VCa_h)                                               : 1

# AHP K+ Current   
    dq/dt = (alphaq - q*(alphaq + 0.001)) /ms                                    : 1
    alphaq = (int(0.00002*Ca<=0.01)*0.00002*Ca + int(0.00002*Ca>0.01)*0.01)      : 1
    I_KAHP = gKAHP_h*q*(Vd-VK_h)                                                 : 1
    
# Ca+-dependent K+ Current    
    dc/dt = (alphac - c*(alphac + betac)) /ms                                    : 1
    alphac = int(Vd <= -10.)*(exp((Vd+50)/11.) - exp((Vd+53.5)/27.)) /18.975\
              + int(Vd > -10.)*2.*exp((-53.5-Vd)/27.)                            : 1
    betac = int(Vd <= -10.)*(2.*exp((-53.5-Vd)/27.)-alphac)                      : 1
    I_KCa = gKCa_h*c*(int(Ca/250.<=1.))*(Vd-VK_h)*Ca/250. \
              + (int(Ca/250.>1.))*gKCa_h*c*(Vd-VK_h)                             : 1

 
# Dendro-Somatic Current    
    I_sd = gsd*(Vs-Vd)                                                           : 1
    
# Intrinsic Parameters
    gCa_h                                                                        : 1
    VL                                                                           : 1
    gL                                                                           : 1
    gsd                                                                          : 1
    Is                                                                           : 1
    Id                                                                           : 1
    I_AMPA_PP                                                                    : 1
    I_GABA_IP                                                                    : 1
'''
#-------------------------------------------------------------------------------

# INTERNEURON Equations---------------------------------------------------------
eqs_IN = '''
    dVi/dt = ((-( I_Li + I_Nai + I_Ki ) - I_AMPA_PI - I_GABA_II + Ii) / Ci) /ms  : 1 
    
# Sodium Current
    mi = alphami / (alphami + betami)                                            : 1
    # alphami = 0.5*(Vi+35.)/(1.-exp(-(Vi+35.)/10.))                             : 1
    alphami = 0.5*10 / exprel(-(Vi+35.)/10.)                                     : 1
    betami = 20.*exp(-(Vi+60.)/18.)                                              : 1
    dhi/dt = (alphahi*(1.-hi)-betahi*hi)  /ms                                    : 1
    alphahi = 0.35*exp(-(Vi+58.)/20.)                                            : 1
    betahi = 5./(1.+exp(-(Vi+28.)/10.))                                          : 1
    I_Nai = gNai*(mi**3)*hi*(Vi-VNai)                                            : 1
    
# Potassium Current
    dni/dt = (alphani*(1.-ni)-betani*ni) /ms                                     : 1
    # alphani = 0.05*(Vi+34.) / (1.-exp(-(Vi+34.)/10.))                          : 1
    alphani = 0.05*10 / exprel(-(Vi+34.)/10.)                                    : 1
    betani = 0.625*exp(-(Vi+44.)/80.)                                            : 1
    I_Ki = gKi*(ni**4)*(Vi-VKi)                                                  : 1
    
# Leakage Current    
    I_Li = gLi*(Vi-VLi)                                                          : 1

#VARYING Network Parameters 
    VLi                                                                          : 1
    gLi                                                                          : 1
    Ii                                                                           : 1
    I_AMPA_PI                                                                    : 1
    I_GABA_II                                                                    : 1
'''
#-------------------------------------------------------------------------------

# SYNAPTIC Equations------------------------------------------------------------
SynAMPA_PP = '''   
    dr_AMPAPP/dt = (Alpha_AMPA * C_AMPAPP * (1 - r_AMPAPP) - Beta_AMPA * r_AMPAPP) /ms           : 1 (clock-driven)
    C_AMPAPP = int(timestep(t - lastupdate_AMPAPP, dt) <= timestep(Cdur_AMPA*ms, dt))* Cmax_AMPA : 1
    Isingle_AMPA_PP = gAMPA_PP * r_AMPAPP * (Vd - Erev_AMPA)                                     : 1 
    I_AMPA_PP_post = Isingle_AMPA_PP                                                             : 1 (summed) 
    gAMPA_PP                                                                                     : 1
    lastupdate_AMPAPP                                                                            : second 
    '''
    
SynGABA_IP = '''  
    dr_GABAIP/dt = (Alpha_GABA * C_GABAIP * (1 - r_GABAIP) - Beta_GABA * r_GABAIP) /ms           : 1 (clock-driven)
    C_GABAIP = int(timestep(t - lastupdate_GABAIP, dt) <= timestep(Cdur_GABA*ms, dt))* Cmax_GABA : 1
    Isingle_GABA_IP = gGABA_IP * r_GABAIP * (Vs - Erev_GABA)                                     : 1 
    I_GABA_IP_post = Isingle_GABA_IP                                                             : 1 (summed)    
    gGABA_IP                                                                                     : 1
    lastupdate_GABAIP                                                                            : second    
'''

SynGABA_II ='''   
    dr_GABAII/dt = (Alpha_GABA * C_GABAII * (1 - r_GABAII) - Beta_GABA * r_GABAII) /ms           : 1 (clock-driven)
    C_GABAII = int(timestep(t - lastupdate_GABAII, dt) <= timestep(Cdur_GABA*ms, dt))* Cmax_GABA : 1
    Isingle_GABA_II = gGABA_II * r_GABAII * (Vi - Erev_GABA)                                       : 1  
    I_GABA_II_post = Isingle_GABA_II                                                             : 1 (summed)
    gGABA_II                                                                                     : 1
    lastupdate_GABAII                                                                            : second  
    '''
    
SynAMPA_PI = '''
    dr_AMPAPI/dt = (Alpha_AMPA * C_AMPAPI * (1 - r_AMPAPI) - Beta_AMPA * r_AMPAPI) /ms           : 1 (clock-driven)
    C_AMPAPI = int(timestep(t - lastupdate_AMPAPI, dt) <= timestep(Cdur_AMPA*ms, dt))* Cmax_AMPA : 1
    Isingle_AMPA_PI = gAMPA_PI * r_AMPAPI * (Vi - Erev_AMPA)                                       : 1 
    I_AMPA_PI_post = Isingle_AMPA_PI                                                             : 1 (summed)
    gAMPA_PI                                                                                     : 1
    lastupdate_AMPAPI                                                                            : second  
    '''
#-------------------------------------------------------------------------------
    
# Equations for when there is a presinaptic spike-------------------------------

GABA_IP_fire = '''
    lastupdate_GABAIP = int(timestep(t - lastupdate_GABAIP, dt) > timestep(Deadtime_GABA*ms, dt)) * t \
    + int(timestep(t - lastupdate_GABAIP, dt) < timestep(Deadtime_GABA*ms, dt)) * lastupdate_GABAIP  
'''

GABA_II_fire = '''
    lastupdate_GABAII = int(timestep(t - lastupdate_GABAII, dt) > timestep(Deadtime_GABA*ms, dt)) * t \
    + int(timestep(t - lastupdate_GABAII, dt) < timestep(Deadtime_GABA*ms, dt)) * lastupdate_GABAII    
'''
AMPA_PP_fire = '''
    lastupdate_AMPAPP = int(timestep(t - lastupdate_AMPAPP, dt) > timestep(Deadtime_AMPA*ms, dt)) * t \
    + int(timestep(t - lastupdate_AMPAPP, dt) < timestep(Deadtime_AMPA*ms, dt)) * lastupdate_AMPAPP    
'''

AMPA_PI_fire = '''
    lastupdate_AMPAPI = int(timestep(t - lastupdate_AMPAPI, dt) > timestep(Deadtime_AMPA*ms, dt)) * t \
    + int(timestep(t - lastupdate_AMPAPI, dt) < timestep(Deadtime_AMPA*ms, dt)) * lastupdate_AMPAPI  
'''
#-------------------------------------------------------------------------------