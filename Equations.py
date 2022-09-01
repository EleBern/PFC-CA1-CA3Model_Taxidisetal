from brian2 import *
from Parameters import *

# PFC PYRAMIDAL CELL Equations--------------------------------------------------
# Changed exponentials with exprel for better accuracy when near division for 0
eqs_PYc = '''
    dVs/dt = (-( I_L + I_Na + I_K + I_A + I_KS + I_KNa ) - I_sd/As -I_GABA_IP + Is) / C_c /ms  : 1 
    dVd/dt = (-( I_NaP + I_AR + I_Ca + I_KCa ) + I_sd/Ad -I_AMPA_PP + Id) / C_c /ms            : 1
    
    # dVs/dt = (-( I_L + I_Na + I_K + I_A + I_KS + I_KNa ) - I_sd/As + Is) / C_c /ms           : 1 
    # dVd/dt = (-( I_NaP + I_AR + I_Ca + I_KCa ) + I_sd/Ad + Id) / C_c /ms       : 1
    
    dNa/dt = ( -alpha_Na*(As*I_Na + Ad*I_NaP)  -  R_pump*( Na**3/(Na**3+15.**3)\
    - Na_eq**3/(Na_eq**3+15.**3) ) )/ms                                          : 1

    dCa/dt = (-alpha_Ca*Ad*I_Ca - Ca/tau_Ca) /ms                                 : 1
    
#SOMATIC currents    
# Sodium Current  
    m = alpham / (alpham + betam)                                                : 1
    dh/dt = phi* (alphah*(1.-h)-betah*h) /ms                                     : 1
    # alpham = 0.1*(Vs+33.) / (1.-exp(-(Vs+33.)/10.))                            : 1
    alpham = 0.1*10 / exprel(-(Vs+33.)/10.)                                      : 1
    betam = 4.*exp(-(Vs+53.7)/12.)                                               : 1
    alphah = 0.07*exp(-(Vs+50.)/10.)                                             : 1
    betah = 1./(1.+exp(-(Vs+20.)/10.))                                           : 1
    I_Na = gNa_c*(m**3)*h*(Vs-VNa_c)                                             : 1

# Potassium Current   
    dn/dt = phi* (alphan*(1.-n)-betan*n) /ms                                     : 1
    # alphan = 0.01*(Vs+34.) / (1.-exp(-(Vs+34.)/10.))                           : 1
    alphan = 0.01*10 / exprel(-(Vs+34.)/10.)                                     : 1
    betan = 0.125*exp(-(Vs+44.)/25.)                                             : 1
    I_K = gK_c*(n**4)*(Vs-VK_c)                                                  : 1
    
# Leakage Current    
    I_L = gL*(Vs-VL)                                                             : 1

# Fast A-type K+ Current    
    mA_inf = 1./(1.+exp(-(Vs+50.)/20.))                                          : 1
    dhA/dt = (hA_inf-hA)/tau_hA    /ms                                           : 1
    hA_inf = 1./(1.+exp((Vs+80.)/6.))                                            : 1
    I_A = gA_c*(mA_inf**3)*hA*(Vs-VK_c)                                          : 1

# Non-inactivating slow K+ Current (Wang 1999a)
    dmKS/dt = (mKS_inf-mKS)/tau_mKS   /ms                                        : 1
    mKS_inf = 1./(1.+exp(-(Vs+34.)/6.5))                                         : 1
    tau_mKS = 8./(exp(-(Vs+55.)/30.) + exp((Vs+55.)/30.))                        : 1 
    I_KS = gKS_c*mKS*(Vs-VK_c)                                                   : 1

# Na+-dependent K+ current  
    w = 0.37/(1.+(38.7/Na)**3.5)                                                 : 1
    I_KNa = gKNa_c*w*(Vs-VK_c)                                                   : 1
    
#DENDRITIC Currents
# Presistent Sodium Current
    mNaP_inf = 1./(1.+exp(-(Vd+55.7)/7.7))                                       : 1
    I_NaP = gNaP_c*(mNaP_inf**3)*(Vd-VNa_c)                                      : 1

# Inward rectifier K+ Current
    hAR_inf = 1./(1.+exp((Vd+75.)/4.))                                           : 1
    I_AR = gAR_c*hAR_inf*(Vd-VK_c)                                               : 1

# High-Threshold Ca+ Current
    mCa_inf = 1./(1.+exp(-(Vd+20.)/9.))                                          : 1
    I_Ca = gCa_c*(mCa_inf**2)*(Vd-VCa_c)                                         : 1

# Ca+-depndent K+ Current    
    I_KCa = gKCa_c*(Ca/(Ca+K_D))*(Vd-VK_c)                                       : 1
 
# Dendro-Somatic Current    
    I_sd = gsd*(Vs-Vd)                                                           : 1

#VARYING Network Paramteres
    VL                                                                           : 1
    gL                                                                           : 1
    gsd                                                                          : 1
    Id                                                                           : 1
    Is                                                                           : 1
    I_AMPA_PP                                                                    : 1
    I_GABA_IP                                                                    : 1
'''
#-------------------------------------------------------------------------------

# INTERNEURON Equations COMPTE--------------------------------------------------
# eqs_IN = '''
#     dVi/dt = ((-( I_Li + I_Nai + I_Ki ) + Ii - I_GABA_II - I_AMPA_PI) / Ci) /ms  : 1
    
# #SOMATIC currents    
# # Sodium Current  
#     mi = alphami / (alphami + betami)                                            : 1
#     dhi/dt = phi* (alphahi*(1.-hi)-betahi*hi) /ms                                : 1
#     alphami = 0.1*(Vi+33.) / (1.-exp(-(Vi+33.)/10.))                             : 1
#     betami = 4.*exp(-(Vi+53.7)/12.)                                              : 1
#     alphahi = 0.07*exp(-(Vi+50.)/10.)                                            : 1
#     betahi = 1./(1.+exp(-(Vi+20.)/10.))                                          : 1
#     I_Nai = gNa_c*(mi**3)*hi*(Vi-VNa_c)                                          : 1

# # Potassium Current   
#     dni/dt = phi* (alphani*(1.-ni)-betani*ni) /ms                                : 1
#     alphani = 0.01*(Vi+34.) / (1.-exp(-(Vi+34.)/10.))                            : 1
#     betani = 0.125*exp(-(Vi+44.)/25.)                                            : 1
#     I_Ki = gK_c*(ni**4)*(Vi-VK_c)                                                : 1
    
# # Leakage Current    
#     I_Li = gL*(Vi-VL)                                                            : 1

# #VARYING Network Parameters 
#     VL                                                                           : 1
#     gL                                                                           : 1
#     Ii                                                                           : 1
# '''


# INTERNEURON Equations TAXIDIS-------------------------------------------------
# Changed exponentials with exprel for better accuracy when near division for 0
eqs_IN = '''
    dVi/dt = ((-( I_Li + I_Nai + I_Ki ) + Ii - I_GABA_II - I_AMPA_PI ) / Ci) /ms : 1 
    # dVi/dt = ((-( I_Li + I_Nai + I_Ki ) + Ii) / Ci) /ms : 1 
    
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
    I_Li = gL*(Vi-VL)                                                            : 1

#VARYING Network Parameters 
    VL                                                                           : 1
    gL                                                                           : 1
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
    Isingle_GABA_II = gGABA_II * r_GABAII * (Vi - Erev_GABA)                                     : 1  
    I_GABA_II_post = Isingle_GABA_II                                                             : 1 (summed)
    gGABA_II                                                                                     : 1
    lastupdate_GABAII                                                                            : second  
    '''
    
SynAMPA_PI = '''
    dr_AMPAPI/dt = (Alpha_AMPA * C_AMPAPI * (1 - r_AMPAPI) - Beta_AMPA * r_AMPAPI) /ms           : 1 (clock-driven)
    C_AMPAPI = int(timestep(t - lastupdate_AMPAPI, dt) <= timestep(Cdur_AMPA*ms, dt))* Cmax_AMPA : 1
    Isingle_AMPA_PI = gAMPA_PI * r_AMPAPI * (Vi - Erev_AMPA)                                     : 1 
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