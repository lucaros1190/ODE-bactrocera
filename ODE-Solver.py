
# Main file of the ODE system solution. Lauch this script to solve the system and
# plot the results.
# Created by Luca Rossini
# e-mail: luca.rossini@unitus.it
# Last update: 17 July 2022


# List of import

from Parameters import *
from ODEbactro import *
import scipy.integrate as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


# Import temperature data and field monitoring

data_exp = pd.read_csv("Experimental_data.csv", sep=";", header=0)
data_exp.columns = ["Day", "Temperature", "Ad_male"]

# Replace NaN values to plot correctly the experimental data

data_exp = data_exp.replace('', np.nan)

# Retrieve the time vector from experimental data

day_ob = data_exp["Day"]
day_obs = day_ob.to_numpy()

# Retrieve the daily temperature array

Temp_data = data_exp['Temperature']
Temp = Temp_data.to_numpy()

# Retrieve monitoring data

yob = data_exp['Ad_male']
yobs = yob.to_numpy() 


# Solve the equation using odeint from scipy - All the stages

y = EqSolver_Full(day_obs, InitCond_Acquired, BriPar_EGG, BriPar_L1, BriPar_L2, BriPar_L3, BriPar_P, BriPar_AD, FertPar_Acquired, MortPar_EGG, MortPar_L1, MortPar_L2, MortPar_L3, MortPar_P, MortPar_AD, SR, Temp)


# Solve the equation using odeint from scipy - Only adult males

y_males = EqSolver(day_obs, InitCond_Acquired, BriPar_EGG, BriPar_L1, BriPar_L2, BriPar_L3, BriPar_P, BriPar_AD, FertPar_Acquired, MortPar_EGG, MortPar_L1, MortPar_L2, MortPar_L3, MortPar_P, MortPar_AD, SR, Temp)


# Plot the ODE system solutions

plt.figure(1)

plt.plot(day_obs, y[:, 0], color="C0", label=f"Egg")
plt.plot(day_obs, y[:, 1], color="C1", label=f"Larva 1")
plt.plot(day_obs, y[:, 2], color="C2", label=f"Larva 2")
plt.plot(day_obs, y[:, 3], color="C3", label=f"Larva 3")
plt.plot(day_obs, y[:, 4], color="C4", label=f"Pupa")
plt.plot(day_obs, y[:, 5], color="C5", label=f"Adult males")
plt.plot(day_obs, y[:, 6], color="C6", label=f"Adult non-mat females")
plt.plot(day_obs, y[:, 7], color="C7", label=f"Adult females")
plt.xlabel('Time (days)')
plt.ylabel('Population density')
plt.legend()


# Plot the rate functions

    # Prepare a linspace temperature array (horizontal axis)

TempForPlot = np.linspace(0, 47, 46)

    # Base temperature-dependent fertility function - Just for plot purposes

def FertFunc_Plot(FertPar, temp):
    
    FP = FertPar[0] * ( ((FertPar[1] + 1) / (np.pi * (FertPar[2] ** (2 * FertPar[1] + 2)))) * ((FertPar[2] ** 2) - ( ((temp - FertPar[4]) ** 2) + (FertPar[3] ** 2)) ) ** FertPar[1])

    return FP 

    # Base temperature-dependent development rate function - Just for plot purposes

def BriFunc_Plot(BriPar, temp):
    
    if np.any(temp >= BriPar[1]) and np.any(temp <= BriPar[2]):
    
        Bri = BriPar[0] * temp * (temp - BriPar[1]) * pow((BriPar[2] - temp), (1 / BriPar[3]))
    else:
        Bri = 0.00001

    return Bri

    # Base temperature-dependent mortality rate funciton - Just for plot purposes

def MortFunc_Plot(MortPar, temp):

    Survival = MortPar[0] * np.exp(1 + ((MortPar[1] - temp) / MortPar[2]) - np.exp((MortPar[1] - temp) / MortPar[2]))

    Mort = 1 - Survival

    return Mort

    # Prepare the data to plot

FertPlot = np.zeros(len(TempForPlot))

for i in range(len(TempForPlot)):
    FertPlot[i] = FertFunc_Plot(FertPar_Acquired, i)

MortPlot_EGG = np.zeros(len(TempForPlot))
MortPlot_L1 = np.zeros(len(TempForPlot))
MortPlot_L2 = np.zeros(len(TempForPlot))
MortPlot_L3 = np.zeros(len(TempForPlot))
MortPlot_P = np.zeros(len(TempForPlot))
MortPlot_AD = np.zeros(len(TempForPlot))
    
for i in range(len(TempForPlot)):
    MortPlot_EGG[i] = MortFunc_Plot(MortPar_EGG, i)
    MortPlot_L1[i] = MortFunc_Plot(MortPar_L1, i)
    MortPlot_L2[i] = MortFunc_Plot(MortPar_L2, i)
    MortPlot_L3[i] = MortFunc_Plot(MortPar_L3, i)
    MortPlot_P[i] = MortFunc_Plot(MortPar_P, i)
    MortPlot_AD[i] = MortFunc_Plot(MortPar_AD, i)


   # Plot the development rate functions - prepare data

DevRate = np.zeros(len(TempForPlot))
    
for i in range(len(TempForPlot)):
    DevRate = BriFunc_Plot(BriPar_EGG, i)


# Plot fertility and mortality rates

fig, FertScatter = plt.subplots()
MortScatter = FertScatter.twinx()

FertScatter.plot(TempForPlot, FertPlot, color="C0", label='Fertility rate')
MortScatter.plot(TempForPlot, MortPlot_EGG, linestyle='-.', color="C1", label='Mortality Egg')
MortScatter.plot(TempForPlot, MortPlot_L1, linestyle='--', color="C2", label='Mortality Larva')
MortScatter.plot(TempForPlot, MortPlot_P, linestyle=':', color="C5", label='Mortality Pupa')
MortScatter.plot(TempForPlot, MortPlot_AD, color="C6", label='Mortality Adults')

FertScatter.set_xlabel('Temperature (°C)')
FertScatter.set_ylabel('Fertility rate (Eggs/female/day)')

MortScatter.set_ylabel('Mortality (Total portion individuals died)')
fig.legend(ncol=2)

    # Plot the development rate functions - plot

plt.figure(4)

plt.plot(TempForPlot, DevRate, color="C0", alpha=0.5, label=f"Egg")
plt.xlabel('Temperature (C)')
plt.ylabel('Development rate (1/day)')
plt.legend()


    # Plot the simulated and experimental adult males

plt.figure(5)

plt.plot(day_obs, y_males, color="C0", label='Simulation')
plt.scatter(day_obs, yobs, color="C1", marker='.', label='Field data')
plt.xlabel('Time (days)')
plt.ylabel('Population density (Number of adult males)')
plt.legend()


# Save the numerical results on a single csv

NumericalResults = pd.DataFrame({'Day': day_obs, 'Egg': y[:, 0], 'L1': y[:, 1], 'L2': y[:, 2], 'L3': y[:, 3], 'Pupa': y[:, 4], 'Ad_male': y[:, 5], 'Ad_NMFemale': y[:, 6], 'Ad_Female': y[:, 7]})

NumericalResults.to_csv('Numerical_results.csv', sep = ';', index = False)


# Make the overall plots

plt.show()



