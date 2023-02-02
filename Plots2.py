# Color plots for matching months and future weather variables
# Atmospheric Innovations Research (AIR) Laboratory, University of Guelph
# Developed by Amir A. Aliabadi and Rachel M. McLeod
# Last updated: 2022-03-16

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
PathFigures = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Figures")
PathResults = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/")
PathEPWFutureCorrected = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/EPW-Future-Corrected")

MatchFutureFile = "Match-CM-EPW-Future.txt"

# Specify the number of past historical, validation, and future years in the dataset
# Specify the number of months in each year
NYearsPast = 20
NYearsValid = 14
NYearsFuture = 80
NMonths = 12
NHoursInDay = 24
NHoursInMonth = 31 * 24
MatchingFileColumns = 6
FirstFutureYear = 2021

CelciusToKelvin = 273.15

# Define starting index for J,F,M,A,M,J,J,A,S,O,N,D in the EPW dataset
# Note: remove February 29 from all the leap years in all the datasets
# Define number of hourly observations in EPW file
MonthIndexCM = [0, 31, 59, 89, 120, 150, 181, 212, 242, 273, 303, 334, 365]
MonthIndexEPW = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
# Note pandas.read_csv() cannot read the first hour of data in the EPW file
MonthIndexEPWPandas = [0, 743, 1415, 2159, 2879, 3623, 4343, 5087, 5831, 6551, 7295, 8015, 8759]
MonthIndex = [744, 672, 744, 720, 744, 720, 744, 744, 720, 744, 720, 744]
EPWNPoints = 8759
DailyNPoints = 365

# Number of rows to skip in the EPW file
NSkipHeader = 8

# Calculations
########################################################################################################################
Years = [i for i in range(FirstFutureYear, FirstFutureYear + NYearsFuture)]

# Loading and analyzing EPW data for the future period
# Go to results directory
os.chdir(PathResults)

# Load matching indices for the future period
MatchingIndicesFuture = numpy.loadtxt(MatchFutureFile)

FutureYear = MatchingIndicesFuture[:,0]
FutureMonth = MatchingIndicesFuture[:,1]
PastYear = MatchingIndicesFuture[:,2]
PastMonth = MatchingIndicesFuture[:,3]

#Define MatchingMonth matrix
FutureMatchingMonth = numpy.zeros((NYearsFuture, NMonths))

#Compute binned data
for i in range(0, NYearsFuture * NMonths):
    FutureYearIndex = int(FutureYear[i])
    FutureMonthIndex = int(FutureMonth[i])
    # Dangerous, because temperature values may be over written in the matrix
    # If there are more than one observation for a given combinations of LatIndex and LonIndex
    # Also depends on the matrix resolution
    FutureMatchingMonth[FutureYearIndex][FutureMonthIndex] = PastMonth[i] + 1

# Go to directory that contains the future EPW file
os.chdir(PathEPWFutureCorrected)
# Remember the file names in order in this directory
EPWFutureFiles = os.listdir(PathEPWFutureCorrected)

# Define matrix for EPW data
EPWFutureTdryb = numpy.zeros((NYearsFuture, EPWNPoints))
EPWFuturePressure = numpy.zeros((NYearsFuture, EPWNPoints))
EPWFutureRad = numpy.zeros((NYearsFuture, EPWNPoints))
EPWFutureWindS = numpy.zeros((NYearsFuture, EPWNPoints))

# Read text files and load data to matrices
i = 0
for file in os.listdir():
    data = pandas.read_csv(file, header=NSkipHeader)
    data = data.values

    # Note pandas.read_csv() cannot read first hour of data in EPW file
    # Consider the second hour data as first hour data
    EPWFutureTdryb[i][0] = data[0, 6]
    EPWFuturePressure[i][0] = data[0, 9]
    # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
    EPWFutureRad[i][0] = data[0, 12] + data[0, 14] + data[0, 15]
    EPWFutureWindS[i][0] = data[0, 21]

    for l in range(1, EPWNPoints):
        EPWFutureTdryb[i][l] = data[l-1, 6]
        EPWFuturePressure[i][l] = data[l-1, 9]
        # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
        EPWFutureRad[i][l] = data[l-1, 12] + data[l-1, 14] + data[l-1, 15]
        EPWFutureWindS[i][l] = data[l-1, 21]
    i = i + 1

print('\nFinished reading historical EPW files.')
print('Last file: ', file)

# Calculate mean and std of weather variables on a monthly basis
EPWMonthlyFutureTdrybMean = numpy.zeros((NYearsFuture, NMonths))
EPWMonthlyFuturePressureMean = numpy.zeros((NYearsFuture, NMonths))
EPWMonthlyFutureRadMean = numpy.zeros((NYearsFuture, NMonths))
EPWMonthlyFutureWindSMean = numpy.zeros((NYearsFuture, NMonths))

EPWMonthlyFutureTdrybStd = numpy.zeros((NYearsFuture, NMonths))
EPWMonthlyFuturePressureStd = numpy.zeros((NYearsFuture, NMonths))
EPWMonthlyFutureRadStd = numpy.zeros((NYearsFuture, NMonths))
EPWMonthlyFutureWindSStd = numpy.zeros((NYearsFuture, NMonths))

for i in range(0, NYearsFuture):
    for j in range(0, NMonths):
        EPWMonthlyFutureTdrybMean[i][j] = numpy.nanmean(EPWFutureTdryb[i][MonthIndexEPW[j]:MonthIndexEPW[j+1]])
        EPWMonthlyFuturePressureMean[i][j] = numpy.nanmean(EPWFuturePressure[i][MonthIndexEPW[j]:MonthIndexEPW[j+1]])
        EPWMonthlyFutureRadMean[i][j] = numpy.nanmean(EPWFutureRad[i][MonthIndexEPW[j]:MonthIndexEPW[j+1]])
        EPWMonthlyFutureWindSMean[i][j] = numpy.nanmean(EPWFutureWindS[i][MonthIndexEPW[j]:MonthIndexEPW[j+1]])

        EPWMonthlyFutureTdrybStd[i][j] = numpy.nanstd(EPWFutureTdryb[i][MonthIndexEPW[j]:MonthIndexEPW[j + 1]])
        EPWMonthlyFuturePressureStd[i][j] = numpy.nanstd(EPWFuturePressure[i][MonthIndexEPW[j]:MonthIndexEPW[j + 1]])
        EPWMonthlyFutureRadStd[i][j] = numpy.nanstd(EPWFutureRad[i][MonthIndexEPW[j]:MonthIndexEPW[j + 1]])
        EPWMonthlyFutureWindSStd[i][j] = numpy.nanstd(EPWFutureWindS[i][MonthIndexEPW[j]:MonthIndexEPW[j + 1]])

TdrybMean = EPWMonthlyFutureTdrybMean
TdrybMeanPStd = EPWMonthlyFutureTdrybMean + EPWMonthlyFutureTdrybStd
TdrybMeanMStd = EPWMonthlyFutureTdrybMean - EPWMonthlyFutureTdrybStd

PressureMean = EPWMonthlyFuturePressureMean
PressureMeanPStd = EPWMonthlyFuturePressureMean + EPWMonthlyFuturePressureStd
PressureMeanMStd = EPWMonthlyFuturePressureMean - EPWMonthlyFuturePressureStd

RadMean = EPWMonthlyFutureRadMean
RadMeanPStd = EPWMonthlyFutureRadMean + EPWMonthlyFutureRadStd
RadMeanMStd = EPWMonthlyFutureRadMean - EPWMonthlyFutureRadStd

WindSMean = EPWMonthlyFutureWindSMean
WindSMeanPStd = EPWMonthlyFutureWindSMean + EPWMonthlyFutureWindSStd
WindSMeanMStd = EPWMonthlyFutureWindSMean - EPWMonthlyFutureWindSStd

# Polyfit analysis to determine trends
for j in range(0, NMonths):
    TdrybMeanFit = numpy.polyfit(Years, EPWMonthlyFutureTdrybMean[:,j], 1)
    print('Month j, TdrybMeanFit: ', j, numpy.round(TdrybMeanFit, 3))

for j in range(0, NMonths):
    PressureMeanFit = numpy.polyfit(Years, EPWMonthlyFuturePressureMean[:, j], 1)
    print('Month j, PressureMeanFit: ', j, numpy.round(PressureMeanFit, 3))

for j in range(0, NMonths):
    RadMeanFit = numpy.polyfit(Years, EPWMonthlyFutureRadMean[:, j], 1)
    print('Month j, RadMeanFit: ', j, numpy.round(RadMeanFit, 3))

for j in range(0, NMonths):
    WindSMeanFit = numpy.polyfit(Years, EPWMonthlyFutureWindSMean[:, j], 1)
    print('Month j, WindSMeanFit: ', j, numpy.round(WindSMeanFit, 3))

# Plots
########################################################################################################################

# Commands to typeset figures using LaTeX formatting
plt.rc('text', usetex = 'True')
plt.rc('font', family = 'serif')
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)

SymbolFormat = ['ko','bs','gd','r<']
LineFormat = ['k-','b-','g-','r-']
ColorFormat = ['k','b','g','r']

# Go to directory of figures
os.chdir(PathFigures)

'''
fig = plt.figure(figsize=(8,6))
for s in range(0, 4):
    plt.plot(Years, TdrybMean[:,3*s], SymbolFormat[s], linewidth = 1)
    plt.plot(Years, TdrybMeanPStd[:,3*s], LineFormat[s], linewidth = 1)
    plt.plot(Years, TdrybMeanMStd[:,3*s], LineFormat[s], linewidth = 1)
    plt.fill_between(Years, TdrybMeanMStd[:,3*s], TdrybMeanPStd[:,3*s], color = ColorFormat[s], alpha=.3)
Jan = mpatches.Patch(color='k', label='January', alpha = 0.5)
Apr = mpatches.Patch(color='b', label='April')
Jul = mpatches.Patch(color='g', label='July')
Oct = mpatches.Patch(color='r', label='October')
plt.legend(handles=[Jan, Apr, Jul, Oct], loc='upper left', fontsize=16)
plt.xlabel('Years', fontsize=16)
plt.ylabel('Dry Bulb Temperature [$^{\circ}$C]', fontsize=16)
fig.show()
'''

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(6,10))
ax1.set_title('January')
ax1.plot(Years, TdrybMean[:,0], SymbolFormat[0], linewidth = 1)
ax1.plot(Years, TdrybMeanPStd[:,0], LineFormat[0], linewidth = 1)
ax1.plot(Years, TdrybMeanMStd[:,0], LineFormat[0], linewidth = 1)
ax1.fill_between(Years, TdrybMeanMStd[:,0], TdrybMeanPStd[:,0], color = ColorFormat[0], alpha=.3)
ax1.set_xticklabels([])
ax1.set_ylim((-20,10))
ax1.set_ylabel('Tdb [$^{\circ}$C]', fontsize=16)
ax2.set_title('April')
ax2.plot(Years, TdrybMean[:,3], SymbolFormat[1], linewidth = 1)
ax2.plot(Years, TdrybMeanPStd[:,3], LineFormat[1], linewidth = 1)
ax2.plot(Years, TdrybMeanMStd[:,3], LineFormat[1], linewidth = 1)
ax2.fill_between(Years, TdrybMeanMStd[:,3], TdrybMeanPStd[:,3], color = ColorFormat[1], alpha=.3)
ax2.set_xticklabels([])
ax2.set_ylim((-10,20))
ax2.set_ylabel('Tdb [$^{\circ}$C]', fontsize=16)
ax3.set_title('July')
ax3.plot(Years, TdrybMean[:,6], SymbolFormat[2], linewidth = 1)
ax3.plot(Years, TdrybMeanPStd[:,6], LineFormat[2], linewidth = 1)
ax3.plot(Years, TdrybMeanMStd[:,6], LineFormat[2], linewidth = 1)
ax3.fill_between(Years, TdrybMeanMStd[:,6], TdrybMeanPStd[:,6], color = ColorFormat[3], alpha=.3)
ax3.set_xticklabels([])
ax3.set_ylim((10,35))
ax3.set_ylabel('Tdb [$^{\circ}$C]', fontsize=16)
ax4.set_title('October')
ax4.plot(Years, TdrybMean[:,9], SymbolFormat[3], linewidth = 1)
ax4.plot(Years, TdrybMeanPStd[:,9], LineFormat[3], linewidth = 1)
ax4.plot(Years, TdrybMeanMStd[:,9], LineFormat[3], linewidth = 1)
ax4.fill_between(Years, TdrybMeanMStd[:,9], TdrybMeanPStd[:,9], color = ColorFormat[3], alpha=.3)
ax4.set_ylim((0,30))
ax4.set_xlabel('Years', fontsize=16)
ax4.set_ylabel('Tdb [$^{\circ}$C]', fontsize=16)
plt.tight_layout()
plt.savefig('FutureTemperatureTimeSeries.pdf', dpi=600)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(6,10))
ax1.set_title('January')
ax1.plot(Years, PressureMean[:,0], SymbolFormat[0], linewidth = 1)
ax1.plot(Years, PressureMeanPStd[:,0], LineFormat[0], linewidth = 1)
ax1.plot(Years, PressureMeanMStd[:,0], LineFormat[0], linewidth = 1)
ax1.fill_between(Years, PressureMeanMStd[:,0], PressureMeanPStd[:,0], color = ColorFormat[0], alpha=.3)
ax1.set_xticklabels([])
ax1.set_ylim((95000,99000))
ax1.set_ylabel('P [Pa]', fontsize=16)
ax2.set_title('April')
ax2.plot(Years, PressureMean[:,3], SymbolFormat[1], linewidth = 1)
ax2.plot(Years, PressureMeanPStd[:,3], LineFormat[1], linewidth = 1)
ax2.plot(Years, PressureMeanMStd[:,3], LineFormat[1], linewidth = 1)
ax2.fill_between(Years, PressureMeanMStd[:,3], PressureMeanPStd[:,3], color = ColorFormat[1], alpha=.3)
ax2.set_xticklabels([])
ax2.set_ylim((95000,99000))
ax2.set_ylabel('P [Pa]', fontsize=16)
ax3.set_title('July')
ax3.plot(Years, PressureMean[:,6], SymbolFormat[2], linewidth = 1)
ax3.plot(Years, PressureMeanPStd[:,6], LineFormat[2], linewidth = 1)
ax3.plot(Years, PressureMeanMStd[:,6], LineFormat[2], linewidth = 1)
ax3.fill_between(Years, PressureMeanMStd[:,6], PressureMeanPStd[:,6], color = ColorFormat[3], alpha=.3)
ax3.set_xticklabels([])
ax3.set_ylim((95000,99000))
ax3.set_ylabel('P [Pa]', fontsize=16)
ax4.set_title('October')
ax4.plot(Years, PressureMean[:,9], SymbolFormat[3], linewidth = 1)
ax4.plot(Years, PressureMeanPStd[:,9], LineFormat[3], linewidth = 1)
ax4.plot(Years, PressureMeanMStd[:,9], LineFormat[3], linewidth = 1)
ax4.fill_between(Years, PressureMeanMStd[:,9], PressureMeanPStd[:,9], color = ColorFormat[3], alpha=.3)
ax4.set_ylim((95000,99000))
ax4.set_xlabel('Years', fontsize=16)
ax4.set_ylabel('P [Pa]', fontsize=16)
plt.tight_layout()
plt.savefig('FuturePressureTimeSeries.pdf', dpi=600)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(6,10))
ax1.set_title('January')
ax1.plot(Years, RadMean[:,0], SymbolFormat[0], linewidth = 1)
ax1.plot(Years, RadMeanPStd[:,0], LineFormat[0], linewidth = 1)
ax1.plot(Years, RadMeanMStd[:,0], LineFormat[0], linewidth = 1)
ax1.fill_between(Years, RadMeanMStd[:,0], RadMeanPStd[:,0], color = ColorFormat[0], alpha=.3)
ax1.set_xticklabels([])
ax1.set_ylim((0,1200))
ax1.set_ylabel('Rad. [W m$^{-2}$]', fontsize=16)
ax2.set_title('April')
ax2.plot(Years, RadMean[:,3], SymbolFormat[1], linewidth = 1)
ax2.plot(Years, RadMeanPStd[:,3], LineFormat[1], linewidth = 1)
ax2.plot(Years, RadMeanMStd[:,3], LineFormat[1], linewidth = 1)
ax2.fill_between(Years, RadMeanMStd[:,3], RadMeanPStd[:,3], color = ColorFormat[1], alpha=.3)
ax2.set_xticklabels([])
ax2.set_ylim((0,1200))
ax2.set_ylabel('Rad. [W m$^{-2}$]', fontsize=16)
ax3.set_title('July')
ax3.plot(Years, RadMean[:,6], SymbolFormat[2], linewidth = 1)
ax3.plot(Years, RadMeanPStd[:,6], LineFormat[2], linewidth = 1)
ax3.plot(Years, RadMeanMStd[:,6], LineFormat[2], linewidth = 1)
ax3.fill_between(Years, RadMeanMStd[:,6], RadMeanPStd[:,6], color = ColorFormat[3], alpha=.3)
ax3.set_xticklabels([])
ax3.set_ylim((0,1200))
ax3.set_ylabel('Rad. [W m$^{-2}$]', fontsize=16)
ax4.set_title('October')
ax4.plot(Years, RadMean[:,9], SymbolFormat[3], linewidth = 1)
ax4.plot(Years, RadMeanPStd[:,9], LineFormat[3], linewidth = 1)
ax4.plot(Years, RadMeanMStd[:,9], LineFormat[3], linewidth = 1)
ax4.fill_between(Years, RadMeanMStd[:,9], RadMeanPStd[:,9], color = ColorFormat[3], alpha=.3)
ax4.set_ylim((0,1200))
ax4.set_xlabel('Years', fontsize=16)
ax4.set_ylabel('Rad. [W m$^{-2}$]', fontsize=16)
plt.tight_layout()
plt.savefig('FutureRadTimeSeries.pdf', dpi=600)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(6,10))
ax1.set_title('January')
ax1.plot(Years, WindSMean[:,0], SymbolFormat[0], linewidth = 1)
ax1.plot(Years, WindSMeanPStd[:,0], LineFormat[0], linewidth = 1)
ax1.plot(Years, WindSMeanMStd[:,0], LineFormat[0], linewidth = 1)
ax1.fill_between(Years, WindSMeanMStd[:,0], WindSMeanPStd[:,0], color = ColorFormat[0], alpha=.3)
ax1.set_xticklabels([])
ax1.set_ylim((0.5,5.5))
ax1.set_ylabel('S [m s$^{-1}$]', fontsize=16)
ax2.set_title('April')
ax2.plot(Years, WindSMean[:,3], SymbolFormat[1], linewidth = 1)
ax2.plot(Years, WindSMeanPStd[:,3], LineFormat[1], linewidth = 1)
ax2.plot(Years, WindSMeanMStd[:,3], LineFormat[1], linewidth = 1)
ax2.fill_between(Years, WindSMeanMStd[:,3], WindSMeanPStd[:,3], color = ColorFormat[1], alpha=.3)
ax2.set_xticklabels([])
ax2.set_ylim((0.5,5.5))
ax2.set_ylabel('S [m s$^{-1}$]', fontsize=16)
ax3.set_title('July')
ax3.plot(Years, WindSMean[:,6], SymbolFormat[2], linewidth = 1)
ax3.plot(Years, WindSMeanPStd[:,6], LineFormat[2], linewidth = 1)
ax3.plot(Years, WindSMeanMStd[:,6], LineFormat[2], linewidth = 1)
ax3.fill_between(Years, WindSMeanMStd[:,6], WindSMeanPStd[:,6], color = ColorFormat[3], alpha=.3)
ax3.set_xticklabels([])
ax3.set_ylim((0.5,5.5))
ax3.set_ylabel('S [m s$^{-1}$]', fontsize=16)
ax4.set_title('October')
ax4.plot(Years, WindSMean[:,9], SymbolFormat[3], linewidth = 1)
ax4.plot(Years, WindSMeanPStd[:,9], LineFormat[3], linewidth = 1)
ax4.plot(Years, WindSMeanMStd[:,9], LineFormat[3], linewidth = 1)
ax4.fill_between(Years, WindSMeanMStd[:,9], WindSMeanPStd[:,9], color = ColorFormat[3], alpha=.3)
ax4.set_ylim((0.5,5.5))
ax4.set_xlabel('Years', fontsize=16)
ax4.set_ylabel('S [m s$^{-1}$]', fontsize=16)
plt.tight_layout()
plt.savefig('FutureWindSTimeSeries.pdf', dpi=600)

FutureMonthAxis = numpy.linspace(1, NMonths + 1, NMonths + 1)
FutureYearAxis = numpy.linspace(FirstFutureYear, FirstFutureYear + NYearsFuture - 1, NYearsFuture)
MonthAxis, YearAxis = numpy.meshgrid(FutureMonthAxis, FutureYearAxis)

fig = plt.figure(figsize=(8,6))
Monthpcolor = plt.pcolor(MonthAxis, YearAxis, FutureMatchingMonth,vmin=1, vmax=12)
cbar = plt.colorbar(Monthpcolor)
cbar.set_label('Matching Month', labelpad=-20, y=1.1, rotation=0, fontsize=16)
ax = plt.axes()
plt.xlabel('Future Month', fontsize=16)
plt.ylabel('Future Year', fontsize=16)
ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
plt.savefig('FutureMatchingMonths.pdf', dpi=600)

plt.show()