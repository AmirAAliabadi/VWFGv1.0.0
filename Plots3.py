# Plot Climate Model (CM) bias for the validation period before and after bias correction
# Developed by Amir A. Aliabadi
# Last updated: 2022-08-09

import numpy
import pandas
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Constants of simulation
########################################################################################################################

# Historical period: 1980-1999
# Validation period: 2007-2020
# Future period: 2021-2100

# Define directories and file names
PathResults = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results")
PathFigures = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Figures")

BiasCMValidNoCorrectionFile = "Bias-CM-EPW-Validation-No-Correction.txt"
BiasCMValidCorrectionFile = "Bias-CM-EPW-Validation-Correction.txt"

# Specify the number of months in each year
NMonths = 12

# Define quantile steps for each month equal to the number of day in month of a non-leap year
QuantileStepsMax = 31
QuantileSteps = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# Calculations
########################################################################################################################

# Read non-bias-corrected climate data for the validation period

# Go to results directory
os.chdir(PathResults)

BiasValidNoCorrectionTdryb = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidNoCorrectionPressure = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidNoCorrectionRad = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidNoCorrectionWindS = numpy.zeros((NMonths, QuantileStepsMax))

DataNoCorrection = numpy.loadtxt(BiasCMValidNoCorrectionFile)
NDataNoCorrection = len(DataNoCorrection)

for i in range(0, NDataNoCorrection):
    BiasValidNoCorrectionTdryb[int(DataNoCorrection[i][0])][int(DataNoCorrection[i][1])] = DataNoCorrection[i][2]
    BiasValidNoCorrectionPressure[int(DataNoCorrection[i][0])][int(DataNoCorrection[i][1])] = DataNoCorrection[i][3]
    BiasValidNoCorrectionRad[int(DataNoCorrection[i][0])][int(DataNoCorrection[i][1])] = DataNoCorrection[i][4]
    BiasValidNoCorrectionWindS[int(DataNoCorrection[i][0])][int(DataNoCorrection[i][1])] = DataNoCorrection[i][5]

BiasValidCorrectionTdryb = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidCorrectionPressure = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidCorrectionRad = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidCorrectionWindS = numpy.zeros((NMonths, QuantileStepsMax))

DataCorrection = numpy.loadtxt(BiasCMValidCorrectionFile)
NDataCorrection = len(DataCorrection)

for i in range(0, NDataCorrection):
    BiasValidCorrectionTdryb[int(DataCorrection[i][0])][int(DataCorrection[i][1])] = DataCorrection[i][2]
    BiasValidCorrectionPressure[int(DataCorrection[i][0])][int(DataCorrection[i][1])] = DataCorrection[i][3]
    BiasValidCorrectionRad[int(DataCorrection[i][0])][int(DataCorrection[i][1])] = DataCorrection[i][4]
    BiasValidCorrectionWindS[int(DataCorrection[i][0])][int(DataCorrection[i][1])] = DataCorrection[i][5]

CDFDaily = numpy.zeros((NMonths, QuantileStepsMax))
for i in range(0, NMonths):
    for j in range(1, QuantileSteps[i] + 1):
        CDFDaily[i][j-1] = j/QuantileSteps[i]

# Plots
########################################################################################################################

# Commands to typeset figures using LaTeX formatting
plt.rc('text', usetex = 'True')
plt.rc('font', family = 'serif')
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)

# Go to directory of figures
os.chdir(PathFigures)

# Plot quantile-quantile data

fig = plt.figure(figsize=(8,6))
plt.plot(BiasValidNoCorrectionTdryb[0][0:QuantileSteps[0]], CDFDaily[0][0:QuantileSteps[0]], 'k-', linewidth=2, label='Jan.')
plt.plot(BiasValidNoCorrectionTdryb[1][0:QuantileSteps[1]], CDFDaily[1][0:QuantileSteps[1]], 'k:', linewidth=2, label='Feb.')
plt.plot(BiasValidNoCorrectionTdryb[2][0:QuantileSteps[2]], CDFDaily[2][0:QuantileSteps[2]], 'k-.', linewidth=2, label='Mar.')
plt.plot(BiasValidNoCorrectionTdryb[3][0:QuantileSteps[3]], CDFDaily[3][0:QuantileSteps[3]], 'g-', linewidth=2, label='Apr.')
plt.plot(BiasValidNoCorrectionTdryb[4][0:QuantileSteps[4]], CDFDaily[4][0:QuantileSteps[4]], 'g:', linewidth=2, label='May.')
plt.plot(BiasValidNoCorrectionTdryb[5][0:QuantileSteps[5]], CDFDaily[5][0:QuantileSteps[5]], 'g-.', linewidth=2, label='Jun.')
plt.plot(BiasValidNoCorrectionTdryb[6][0:QuantileSteps[6]], CDFDaily[6][0:QuantileSteps[6]], 'r-', linewidth=2, label='Jul.')
plt.plot(BiasValidNoCorrectionTdryb[7][0:QuantileSteps[7]], CDFDaily[7][0:QuantileSteps[7]], 'r:', linewidth=2, label='Aug.')
plt.plot(BiasValidNoCorrectionTdryb[8][0:QuantileSteps[8]], CDFDaily[8][0:QuantileSteps[8]], 'r-.', linewidth=2, label='Sep.')
plt.plot(BiasValidNoCorrectionTdryb[9][0:QuantileSteps[9]], CDFDaily[9][0:QuantileSteps[9]], 'b-', linewidth=2, label='Oct.')
plt.plot(BiasValidNoCorrectionTdryb[10][0:QuantileSteps[10]], CDFDaily[10][0:QuantileSteps[10]], 'b:', linewidth=2, label='Nov.')
plt.plot(BiasValidNoCorrectionTdryb[11][0:QuantileSteps[11]], CDFDaily[11][0:QuantileSteps[11]], 'b-.', linewidth=2, label='Dec.')
plt.xlabel('Dry Bulb Temperature Bias [$^{\circ}$C]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
plt.axvline(linewidth=2, color='grey')
plt.legend(fontsize=19)
plt.xlim((-10,10))
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQTdrybBiasCMNoCorrection.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
plt.plot(BiasValidCorrectionTdryb[0][0:QuantileSteps[0]], CDFDaily[0][0:QuantileSteps[0]], 'k-', linewidth=2, label='Jan.')
plt.plot(BiasValidCorrectionTdryb[1][0:QuantileSteps[1]], CDFDaily[1][0:QuantileSteps[1]], 'k:', linewidth=2, label='Feb.')
plt.plot(BiasValidCorrectionTdryb[2][0:QuantileSteps[2]], CDFDaily[2][0:QuantileSteps[2]], 'k-.', linewidth=2, label='Mar.')
plt.plot(BiasValidCorrectionTdryb[3][0:QuantileSteps[3]], CDFDaily[3][0:QuantileSteps[3]], 'g-', linewidth=2, label='Apr.')
plt.plot(BiasValidCorrectionTdryb[4][0:QuantileSteps[4]], CDFDaily[4][0:QuantileSteps[4]], 'g:', linewidth=2, label='May.')
plt.plot(BiasValidCorrectionTdryb[5][0:QuantileSteps[5]], CDFDaily[5][0:QuantileSteps[5]], 'g-.', linewidth=2, label='Jun.')
plt.plot(BiasValidCorrectionTdryb[6][0:QuantileSteps[6]], CDFDaily[6][0:QuantileSteps[6]], 'r-', linewidth=2, label='Jul.')
plt.plot(BiasValidCorrectionTdryb[7][0:QuantileSteps[7]], CDFDaily[7][0:QuantileSteps[7]], 'r:', linewidth=2, label='Aug.')
plt.plot(BiasValidCorrectionTdryb[8][0:QuantileSteps[8]], CDFDaily[8][0:QuantileSteps[8]], 'r-.', linewidth=2, label='Sep.')
plt.plot(BiasValidCorrectionTdryb[9][0:QuantileSteps[9]], CDFDaily[9][0:QuantileSteps[9]], 'b-', linewidth=2, label='Oct.')
plt.plot(BiasValidCorrectionTdryb[10][0:QuantileSteps[10]], CDFDaily[10][0:QuantileSteps[10]], 'b:', linewidth=2, label='Nov.')
plt.plot(BiasValidCorrectionTdryb[11][0:QuantileSteps[11]], CDFDaily[11][0:QuantileSteps[11]], 'b-.', linewidth=2, label='Dec.')
plt.xlabel('Dry Bulb Temperature Bias [$^{\circ}$C]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
plt.axvline(linewidth=2, color='grey')
plt.legend(fontsize=19)
plt.xlim((-10,10))
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQTdrybBiasCMCorrection.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
plt.plot(BiasValidNoCorrectionPressure[0][0:QuantileSteps[0]], CDFDaily[0][0:QuantileSteps[0]], 'k-', linewidth=2, label='Jan.')
plt.plot(BiasValidNoCorrectionPressure[1][0:QuantileSteps[1]], CDFDaily[1][0:QuantileSteps[1]], 'k:', linewidth=2, label='Feb.')
plt.plot(BiasValidNoCorrectionPressure[2][0:QuantileSteps[2]], CDFDaily[2][0:QuantileSteps[2]], 'k-.', linewidth=2, label='Mar.')
plt.plot(BiasValidNoCorrectionPressure[3][0:QuantileSteps[3]], CDFDaily[3][0:QuantileSteps[3]], 'g-', linewidth=2, label='Apr.')
plt.plot(BiasValidNoCorrectionPressure[4][0:QuantileSteps[4]], CDFDaily[4][0:QuantileSteps[4]], 'g:', linewidth=2, label='May.')
plt.plot(BiasValidNoCorrectionPressure[5][0:QuantileSteps[5]], CDFDaily[5][0:QuantileSteps[5]], 'g-.', linewidth=2, label='Jun.')
plt.plot(BiasValidNoCorrectionPressure[6][0:QuantileSteps[6]], CDFDaily[6][0:QuantileSteps[6]], 'r-', linewidth=2, label='Jul.')
plt.plot(BiasValidNoCorrectionPressure[7][0:QuantileSteps[7]], CDFDaily[7][0:QuantileSteps[7]], 'r:', linewidth=2, label='Aug.')
plt.plot(BiasValidNoCorrectionPressure[8][0:QuantileSteps[8]], CDFDaily[8][0:QuantileSteps[8]], 'r-.', linewidth=2, label='Sep.')
plt.plot(BiasValidNoCorrectionPressure[9][0:QuantileSteps[9]], CDFDaily[9][0:QuantileSteps[9]], 'b-', linewidth=2, label='Oct.')
plt.plot(BiasValidNoCorrectionPressure[10][0:QuantileSteps[10]], CDFDaily[10][0:QuantileSteps[10]], 'b:', linewidth=2, label='Nov.')
plt.plot(BiasValidNoCorrectionPressure[11][0:QuantileSteps[11]], CDFDaily[11][0:QuantileSteps[11]], 'b-.', linewidth=2, label='Dec.')
plt.xlabel('Pressure Bias [Pa]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
plt.axvline(linewidth=2, color='grey')
plt.legend(fontsize=19)
plt.xlim((-2000,2000))
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQPressureBiasCMNoCorrection.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
plt.plot(BiasValidCorrectionPressure[0][0:QuantileSteps[0]], CDFDaily[0][0:QuantileSteps[0]], 'k-', linewidth=2, label='Jan.')
plt.plot(BiasValidCorrectionPressure[1][0:QuantileSteps[1]], CDFDaily[1][0:QuantileSteps[1]], 'k:', linewidth=2, label='Feb.')
plt.plot(BiasValidCorrectionPressure[2][0:QuantileSteps[2]], CDFDaily[2][0:QuantileSteps[2]], 'k-.', linewidth=2, label='Mar.')
plt.plot(BiasValidCorrectionPressure[3][0:QuantileSteps[3]], CDFDaily[3][0:QuantileSteps[3]], 'g-', linewidth=2, label='Apr.')
plt.plot(BiasValidCorrectionPressure[4][0:QuantileSteps[4]], CDFDaily[4][0:QuantileSteps[4]], 'g:', linewidth=2, label='May.')
plt.plot(BiasValidCorrectionPressure[5][0:QuantileSteps[5]], CDFDaily[5][0:QuantileSteps[5]], 'g-.', linewidth=2, label='Jun.')
plt.plot(BiasValidCorrectionPressure[6][0:QuantileSteps[6]], CDFDaily[6][0:QuantileSteps[6]], 'r-', linewidth=2, label='Jul.')
plt.plot(BiasValidCorrectionPressure[7][0:QuantileSteps[7]], CDFDaily[7][0:QuantileSteps[7]], 'r:', linewidth=2, label='Aug.')
plt.plot(BiasValidCorrectionPressure[8][0:QuantileSteps[8]], CDFDaily[8][0:QuantileSteps[8]], 'r-.', linewidth=2, label='Sep.')
plt.plot(BiasValidCorrectionPressure[9][0:QuantileSteps[9]], CDFDaily[9][0:QuantileSteps[9]], 'b-', linewidth=2, label='Oct.')
plt.plot(BiasValidCorrectionPressure[10][0:QuantileSteps[10]], CDFDaily[10][0:QuantileSteps[10]], 'b:', linewidth=2, label='Nov.')
plt.plot(BiasValidCorrectionPressure[11][0:QuantileSteps[11]], CDFDaily[11][0:QuantileSteps[11]], 'b-.', linewidth=2, label='Dec.')
plt.xlabel('Pressure Bias [Pa]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
plt.axvline(linewidth=2, color='grey')
plt.legend(fontsize=19)
plt.xlim((-2000,2000))
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQPressureBiasCMCorrection.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
plt.plot(BiasValidNoCorrectionRad[0][0:QuantileSteps[0]], CDFDaily[0][0:QuantileSteps[0]], 'k-', linewidth=2, label='Jan.')
plt.plot(BiasValidNoCorrectionRad[1][0:QuantileSteps[1]], CDFDaily[1][0:QuantileSteps[1]], 'k:', linewidth=2, label='Feb.')
plt.plot(BiasValidNoCorrectionRad[2][0:QuantileSteps[2]], CDFDaily[2][0:QuantileSteps[2]], 'k-.', linewidth=2, label='Mar.')
plt.plot(BiasValidNoCorrectionRad[3][0:QuantileSteps[3]], CDFDaily[3][0:QuantileSteps[3]], 'g-', linewidth=2, label='Apr.')
plt.plot(BiasValidNoCorrectionRad[4][0:QuantileSteps[4]], CDFDaily[4][0:QuantileSteps[4]], 'g:', linewidth=2, label='May.')
plt.plot(BiasValidNoCorrectionRad[5][0:QuantileSteps[5]], CDFDaily[5][0:QuantileSteps[5]], 'g-.', linewidth=2, label='Jun.')
plt.plot(BiasValidNoCorrectionRad[6][0:QuantileSteps[6]], CDFDaily[6][0:QuantileSteps[6]], 'r-', linewidth=2, label='Jul.')
plt.plot(BiasValidNoCorrectionRad[7][0:QuantileSteps[7]], CDFDaily[7][0:QuantileSteps[7]], 'r:', linewidth=2, label='Aug.')
plt.plot(BiasValidNoCorrectionRad[8][0:QuantileSteps[8]], CDFDaily[8][0:QuantileSteps[8]], 'r-.', linewidth=2, label='Sep.')
plt.plot(BiasValidNoCorrectionRad[9][0:QuantileSteps[9]], CDFDaily[9][0:QuantileSteps[9]], 'b-', linewidth=2, label='Oct.')
plt.plot(BiasValidNoCorrectionRad[10][0:QuantileSteps[10]], CDFDaily[10][0:QuantileSteps[10]], 'b:', linewidth=2, label='Nov.')
plt.plot(BiasValidNoCorrectionRad[11][0:QuantileSteps[11]], CDFDaily[11][0:QuantileSteps[11]], 'b-.', linewidth=2, label='Dec.')
plt.xlabel('Global Horiz. Rad. Bias [W m$^{-2}$]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
plt.axvline(linewidth=2, color='grey')
plt.legend(fontsize=19)
plt.xlim((-100,100))
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQRadBiasCMNoCorrection.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
plt.plot(BiasValidCorrectionRad[0][0:QuantileSteps[0]], CDFDaily[0][0:QuantileSteps[0]], 'k-', linewidth=2, label='Jan.')
plt.plot(BiasValidCorrectionRad[1][0:QuantileSteps[1]], CDFDaily[1][0:QuantileSteps[1]], 'k:', linewidth=2, label='Feb.')
plt.plot(BiasValidCorrectionRad[2][0:QuantileSteps[2]], CDFDaily[2][0:QuantileSteps[2]], 'k-.', linewidth=2, label='Mar.')
plt.plot(BiasValidCorrectionRad[3][0:QuantileSteps[3]], CDFDaily[3][0:QuantileSteps[3]], 'g-', linewidth=2, label='Apr.')
plt.plot(BiasValidCorrectionRad[4][0:QuantileSteps[4]], CDFDaily[4][0:QuantileSteps[4]], 'g:', linewidth=2, label='May.')
plt.plot(BiasValidCorrectionRad[5][0:QuantileSteps[5]], CDFDaily[5][0:QuantileSteps[5]], 'g-.', linewidth=2, label='Jun.')
plt.plot(BiasValidCorrectionRad[6][0:QuantileSteps[6]], CDFDaily[6][0:QuantileSteps[6]], 'r-', linewidth=2, label='Jul.')
plt.plot(BiasValidCorrectionRad[7][0:QuantileSteps[7]], CDFDaily[7][0:QuantileSteps[7]], 'r:', linewidth=2, label='Aug.')
plt.plot(BiasValidCorrectionRad[8][0:QuantileSteps[8]], CDFDaily[8][0:QuantileSteps[8]], 'r-.', linewidth=2, label='Sep.')
plt.plot(BiasValidCorrectionRad[9][0:QuantileSteps[9]], CDFDaily[9][0:QuantileSteps[9]], 'b-', linewidth=2, label='Oct.')
plt.plot(BiasValidCorrectionRad[10][0:QuantileSteps[10]], CDFDaily[10][0:QuantileSteps[10]], 'b:', linewidth=2, label='Nov.')
plt.plot(BiasValidCorrectionRad[11][0:QuantileSteps[11]], CDFDaily[11][0:QuantileSteps[11]], 'b-.', linewidth=2, label='Dec.')
plt.xlabel('Global Horiz. Rad. Bias [W m$^{-2}$]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
plt.axvline(linewidth=2, color='grey')
plt.legend(fontsize=19)
plt.xlim((-100,100))
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQRadBiasCMCorrection.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
plt.plot(BiasValidNoCorrectionWindS[0][0:QuantileSteps[0]], CDFDaily[0][0:QuantileSteps[0]], 'k-', linewidth=2, label='Jan.')
plt.plot(BiasValidNoCorrectionWindS[1][0:QuantileSteps[1]], CDFDaily[1][0:QuantileSteps[1]], 'k:', linewidth=2, label='Feb.')
plt.plot(BiasValidNoCorrectionWindS[2][0:QuantileSteps[2]], CDFDaily[2][0:QuantileSteps[2]], 'k-.', linewidth=2, label='Mar.')
plt.plot(BiasValidNoCorrectionWindS[3][0:QuantileSteps[3]], CDFDaily[3][0:QuantileSteps[3]], 'g-', linewidth=2, label='Apr.')
plt.plot(BiasValidNoCorrectionWindS[4][0:QuantileSteps[4]], CDFDaily[4][0:QuantileSteps[4]], 'g:', linewidth=2, label='May.')
plt.plot(BiasValidNoCorrectionWindS[5][0:QuantileSteps[5]], CDFDaily[5][0:QuantileSteps[5]], 'g-.', linewidth=2, label='Jun.')
plt.plot(BiasValidNoCorrectionWindS[6][0:QuantileSteps[6]], CDFDaily[6][0:QuantileSteps[6]], 'r-', linewidth=2, label='Jul.')
plt.plot(BiasValidNoCorrectionWindS[7][0:QuantileSteps[7]], CDFDaily[7][0:QuantileSteps[7]], 'r:', linewidth=2, label='Aug.')
plt.plot(BiasValidNoCorrectionWindS[8][0:QuantileSteps[8]], CDFDaily[8][0:QuantileSteps[8]], 'r-.', linewidth=2, label='Sep.')
plt.plot(BiasValidNoCorrectionWindS[9][0:QuantileSteps[9]], CDFDaily[9][0:QuantileSteps[9]], 'b-', linewidth=2, label='Oct.')
plt.plot(BiasValidNoCorrectionWindS[10][0:QuantileSteps[10]], CDFDaily[10][0:QuantileSteps[10]], 'b:', linewidth=2, label='Nov.')
plt.plot(BiasValidNoCorrectionWindS[11][0:QuantileSteps[11]], CDFDaily[11][0:QuantileSteps[11]], 'b-.', linewidth=2, label='Dec.')
plt.xlabel('Wind Speed Bias [m s$^{-1}$]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
plt.axvline(linewidth=2, color='grey')
plt.legend(fontsize=19)
plt.xlim((-2,2))
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQWindSBiasCMNoCorrection.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
plt.plot(BiasValidCorrectionWindS[0][0:QuantileSteps[0]], CDFDaily[0][0:QuantileSteps[0]], 'k-', linewidth=2, label='Jan.')
plt.plot(BiasValidCorrectionWindS[1][0:QuantileSteps[1]], CDFDaily[1][0:QuantileSteps[1]], 'k:', linewidth=2, label='Feb.')
plt.plot(BiasValidCorrectionWindS[2][0:QuantileSteps[2]], CDFDaily[2][0:QuantileSteps[2]], 'k-.', linewidth=2, label='Mar.')
plt.plot(BiasValidCorrectionWindS[3][0:QuantileSteps[3]], CDFDaily[3][0:QuantileSteps[3]], 'g-', linewidth=2, label='Apr.')
plt.plot(BiasValidCorrectionWindS[4][0:QuantileSteps[4]], CDFDaily[4][0:QuantileSteps[4]], 'g:', linewidth=2, label='May.')
plt.plot(BiasValidCorrectionWindS[5][0:QuantileSteps[5]], CDFDaily[5][0:QuantileSteps[5]], 'g-.', linewidth=2, label='Jun.')
plt.plot(BiasValidCorrectionWindS[6][0:QuantileSteps[6]], CDFDaily[6][0:QuantileSteps[6]], 'r-', linewidth=2, label='Jul.')
plt.plot(BiasValidCorrectionWindS[7][0:QuantileSteps[7]], CDFDaily[7][0:QuantileSteps[7]], 'r:', linewidth=2, label='Aug.')
plt.plot(BiasValidCorrectionWindS[8][0:QuantileSteps[8]], CDFDaily[8][0:QuantileSteps[8]], 'r-.', linewidth=2, label='Sep.')
plt.plot(BiasValidCorrectionWindS[9][0:QuantileSteps[9]], CDFDaily[9][0:QuantileSteps[9]], 'b-', linewidth=2, label='Oct.')
plt.plot(BiasValidCorrectionWindS[10][0:QuantileSteps[10]], CDFDaily[10][0:QuantileSteps[10]], 'b:', linewidth=2, label='Nov.')
plt.plot(BiasValidCorrectionWindS[11][0:QuantileSteps[11]], CDFDaily[11][0:QuantileSteps[11]], 'b-.', linewidth=2, label='Dec.')
plt.xlabel('Wind Speed Bias [m s$^{-1}$]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
plt.axvline(linewidth=2, color='grey')
plt.legend(fontsize=19)
plt.xlim((-2,2))
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQWindSBiasCMCorrection.pdf', dpi=600)
fig.show()

plt.show()