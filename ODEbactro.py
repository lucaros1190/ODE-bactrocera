# Python script to solve the ODE model applied to Bactrocera oleae
# This file contains all the functions needed
# Created by Luca Rossini
# e-mail: luca.rossini@unitus.it
# Last update: 29 January 2022


# List of import

import pandas as pd
import numpy as np
from Parameters import *
import scipy.integrate as sc


# Parameter organization

    # Rate functions parameters from the file 'Parameters.py'

BriPar_EGG = a_E, T_L_E, T_M_E, m_E
BriPar_L1 = a_L1, T_L_L1, T_M_L1, m_L1
BriPar_L2 = a_L2, T_L_L2, T_M_L2, m_L2
BriPar_L3 = a_L3, T_L_L3, T_M_L3, m_L3
BriPar_P = a_P, T_L_P, T_M_P, m_P
BriPar_AD = a_AD, T_L_AD, T_M_AD, m_AD

FertPar_Acquired = alpha, gamma, Lambda, delta, tau

MortPar_EGG = k_E, T_MAX_E, rho_E
MortPar_L1 = k_L1, T_MAX_L1, rho_L1
MortPar_L2 = k_L2, T_MAX_L2, rho_L2
MortPar_L3 = k_L3, T_MAX_L3, rho_L3
MortPar_P = k_P, T_MAX_P, rho_P
MortPar_AD = k_AD, T_MAX_AD, rho_AD

    # Initial conditions from the file 'Parameters.py'

InitCond_Acquired = E_0, L1_0, L2_0, L3_0, P_0, Am_0, Anmf_0, Amf_0


# Definition of Briere rate function

def BriFunc(BriPar, DailyTemp, time):
    
    try:
        temp = DailyTemp[int(time)]
    except IndexError:
        print('Warning: Skipping IndexError in BriFunc()! \n')
        temp = DailyTemp[0]
    except KeyError:
        print('Warning: Skipping KeyError in BriFunc()! \n')
        temp = DailyTemp[0]
    
    if np.any(temp >= BriPar[1]) and np.any(temp <= BriPar[2]):
    
        Bri = BriPar[0] * temp * (temp - BriPar[1]) * pow((BriPar[2] - temp), (1 / BriPar[3]))
    else:
        Bri = 0.00001

    return Bri


# Definition of the fertility rate function

def FertFunc(FertPar, DailyTemp, time):

    try:
        temp = DailyTemp[int(time)]
    except IndexError:
        print('Warning: Skipping IndexError in FertFunc()! \n')
        temp = DailyTemp[0]
    except KeyError:
        print('Warning: Skipping KeyError in FertFunc()! \n')
        temp = DailyTemp[0]
    
    FP = FertPar[0] * ( ((FertPar[1] + 1) / (np.pi * (FertPar[2] ** (2 * FertPar[1] + 2)))) * ((FertPar[2] ** 2) - ( ((temp - FertPar[4]) ** 2) + (FertPar[3] ** 2)) ) ** FertPar[1])

    return FP


# Definition of the Mortality rate function as complementary of the survival rate

def MortFunc(MortPar, DailyTemp, time):

    try:
        temp = DailyTemp[int(time)]
    except IndexError:
        print('Warning: Skipping IndexError in MortFunc()! \n')
        temp = DailyTemp[0]

    Survival = MortPar[0] * np.exp(1 + ((MortPar[1] - temp) / MortPar[2]) - np.exp((MortPar[1] - temp) / MortPar[2]))

    Mort = 1 - Survival

    return Mort 


# Definition of the ODE system - odeint function

def Sys_ODE(y, time, BriPar_EGG, BriPar_L1, BriPar_L2, BriPar_L3, BriPar_P, BriPar_AD, FertPar, MortPar_EGG, MortPar_L1, MortPar_L2, MortPar_L3, MortPar_P, MortPar_AD, SexRatio, DailyTemp):

    E = FertFunc(FertPar, DailyTemp, time) * y[7] - BriFunc(BriPar_EGG, DailyTemp, time) * y[0] - (BriFunc(BriPar_EGG, DailyTemp, time) * MortFunc(MortPar_EGG, DailyTemp, time)) * y[0]

    L1 = BriFunc(BriPar_EGG, DailyTemp, time) * y[0] - BriFunc(BriPar_L1, DailyTemp, time) * y[1] - (BriFunc(BriPar_L1, DailyTemp, time) * MortFunc(MortPar_L1, DailyTemp, time)) * y[1]

    L2 = BriFunc(BriPar_L1, DailyTemp, time) * y[1] - BriFunc(BriPar_L2, DailyTemp, time) * y[2] - (BriFunc(BriPar_L2, DailyTemp, time) * MortFunc(MortPar_L2, DailyTemp, time)) * y[2]

    L3 = BriFunc(BriPar_L2, DailyTemp, time) * y[2] - BriFunc(BriPar_L3, DailyTemp, time) * y[3] - (BriFunc(BriPar_L3, DailyTemp, time) * MortFunc(MortPar_L3, DailyTemp, time)) * y[3]

    P = BriFunc(BriPar_L3, DailyTemp, time) * y[3] - BriFunc(BriPar_P, DailyTemp, time) * y[4] - (BriFunc(BriPar_P, DailyTemp, time) * MortFunc(MortPar_P, DailyTemp, time)) * y[4]

    Am = (1 - SexRatio) * BriFunc(BriPar_P, DailyTemp, time) * y[4] - BriFunc(BriPar_AD, DailyTemp, time) * y[5] - (BriFunc(BriPar_AD, DailyTemp, time) * MortFunc(MortPar_AD, DailyTemp, time)) * y[5]

    Anmf = SexRatio * BriFunc(BriPar_P, DailyTemp, time) * y[4] - y[6]

    Amf = y[6] - BriFunc(BriPar_AD, DailyTemp, time) * y[6] - (BriFunc(BriPar_AD, DailyTemp, time) * MortFunc(MortPar_AD, DailyTemp, time)) * y[6] - (BriFunc(BriPar_AD, DailyTemp, time) * MortFunc(MortPar_AD, DailyTemp, time)) * y[7] - BriFunc(BriPar_AD, DailyTemp, time) * y[7]

    Z = np.array([E, L1, L2, L3, P, Am, Anmf, Amf])
    
    return Z


# This function feeds the parameters and solves the equation - It returns the ODE solution of the adult males for comparison

# NOTE: if you want to change the stage to compare, change the output of this function!!

def EqSolver(time_Eq, InitCond_Eq, BriPar_Eq_EGG, BriPar_Eq_L1, BriPar_Eq_L2, BriPar_Eq_L3, BriPar_Eq_P, BriPar_Eq_AD, FertPar_Eq, MortPar_Eq_EGG, MortPar_Eq_L1, MortPar_Eq_L2, MortPar_Eq_L3, MortPar_Eq_P, MortPar_Eq_AD, SR_Eq, Temp_Eq):

    ODE_Sol = sc.odeint(Sys_ODE, y0 = InitCond_Eq, t = time_Eq, args = (BriPar_Eq_EGG, BriPar_Eq_L1, BriPar_Eq_L2, BriPar_Eq_L3, BriPar_Eq_P, BriPar_Eq_AD, FertPar_Eq, MortPar_Eq_EGG, MortPar_Eq_L1, MortPar_Eq_L2, MortPar_Eq_L3, MortPar_Eq_P, MortPar_Eq_AD, SR_Eq, Temp_Eq,))

    return ODE_Sol[:, 5]


# This function feeds the parameters and solves the equation - It returns the ODE solution of all the system

def EqSolver_Full(time_Eq, InitCond_Eq, BriPar_Eq_EGG, BriPar_Eq_L1, BriPar_Eq_L2, BriPar_Eq_L3, BriPar_Eq_P, BriPar_Eq_AD, FertPar_Eq, MortPar_Eq_EGG, MortPar_Eq_L1, MortPar_Eq_L2, MortPar_Eq_L3, MortPar_Eq_P, MortPar_Eq_AD, SR_Eq, Temp_Eq):

    ODE_Sol = sc.odeint(Sys_ODE, y0 = InitCond_Eq, t = time_Eq, args = (BriPar_Eq_EGG, BriPar_Eq_L1, BriPar_Eq_L2, BriPar_Eq_L3, BriPar_Eq_P, BriPar_Eq_AD, FertPar_Eq, MortPar_Eq_EGG, MortPar_Eq_L1, MortPar_Eq_L2, MortPar_Eq_L3, MortPar_Eq_P, MortPar_Eq_AD, SR_Eq, Temp_Eq,))

    return ODE_Sol
 


