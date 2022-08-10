
# List of the parameters in input into the script 'ODEbactro.py'


# Fertility

    # Temperature-dependence

alpha = 7000.0
gamma = 68.444549
Lambda = 50.0
delta = 7.1315087
tau = 25.6458002


# Development rates

    # EGG - Briere 

a_E = 4.60 * pow(10, -5)
T_L_E = 7.23
T_M_E = 32.2
m_E = 2.5

    # L1 - Briere 

a_L1 = 4.60 * pow(10, -5)
T_L_L1 = 7.23
T_M_L1 = 32.2
m_L1 = 2.5

    # L2 - Briere 

a_L2 = 4.60 * pow(10, -5)
T_L_L2 = 7.23
T_M_L2 = 32.2
m_L2 = 2.5

    # L3 - Briere 

a_L3 = 4.60 * pow(10, -5)
T_L_L3 = 7.23
T_M_L3 = 32.2
m_L3 = 2.5

    # P - Briere 

a_P = 4.60 * pow(10, -5)
T_L_P = 7.23
T_M_P = 32.2
m_P = 2.5

    # ADULT - Briere 

a_AD = 4.60 * pow(10, -5)
T_L_AD = 7.23
T_M_AD = 32.2
m_AD = 2.5


# Survival/Mortality rates 

    # EGG - Mortality (1 - survival)

k_E = 0.8660965
T_MAX_E = 25.0
rho_E = -4.5

    # L1 - Mortality (1 - survival)

k_L1 = 0.8515692
T_MAX_L1 = 22.2953273 
rho_L1 = -4.53

    # L2 - Mortality (1 - survival)

k_L2 = 0.8515692
T_MAX_L2 = 22.2953273 
rho_L2 = -4.53

    # L3 - Mortality (1 - survival)

k_L3 = 0.8515692
T_MAX_L3 = 22.2953273 
rho_L3 = -4.53

    # P - Mortality (1 - survival)

k_P = 0.82454
T_MAX_P = 24.6104718
rho_P = -5.9965168

    # ADULT - Mortality (1 - survival)

k_AD = 0.6383781
T_MAX_AD = 24.9275845
rho_AD = -5.9619852


# Sex ratio

SR = 0.5


# Initial conditions for the ODE system

E_0 = 0 #0 <-- This is the stage number into the y[] of Sys_ODE
L1_0 = 0 #1
L2_0 = 0 #2
L3_0 = 0 #3
P_0 = 0 #4
Am_0 = 0 #5
Anmf_0 = 0 #6
Amf_0 = 0.2 #7

