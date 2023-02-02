# Quantile-Quantile plots for all historical records and a sample validation record
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
PathEPWHistorical = ("C:/GoogleDrive/U Guelph/Projects/VWFG/EPW/EPW-Historical")
PathEPWValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/EPW/EPW-Validation")
PathCMValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/CM-Validation-Bias-Corrected")

# Specify the number of past historical, validation, and future years in the dataset
# Specify the number of months in each year
NYearsPast = 20
NYearsValid = 14
NMonths = 12
NHoursInDay = 24

CelciusToKelvin = 273.15

# Define starting index for J,F,M,A,M,J,J,A,S,O,N,D in the EPW dataset
# Note: remove February 29 from all the leap years in all the datasets
# Define number of hourly observations in EPW file
MonthIndex = [0, 31, 59, 89, 120, 150, 181, 212, 242, 273, 303, 334]
EPWNPoints = 8760
DailyNPoints = 365

# Number of rows to skip in the EPW file
NSkipHeader = 8

# Define quantile steps for each month equal to the number of day in month of a non-leap year
QuantileStepsMax = 31
QuantileSteps = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

PlotYearValidIndex = 2
PlotMonthValidIndex = 11

# Define weights for each month and weather variable according to Hosseini et al. (2021) for Montreal
#'''
FSWeightsTdryb = [0.773, 0.912, 0.891, 0.884, 0.869, 0.887, 0.830, 0.868, 0.829, 0.831, 0.905, 0.922]
FSWeightsPressure = [0.086, 0.036, 0.041, 0.043, 0.052, 0.048, 0.064, 0.057, 0.068, 0.060, 0.039, 0.027]
FSWeightsRad = [0.078, 0.027, 0.036, 0.038, 0.043, 0.037, 0.057, 0.045, 0.054, 0.056, 0.031, 0.026]
FSWeightsWindS = [0.063, 0.026, 0.032, 0.036, 0.036, 0.028, 0.048, 0.030, 0.048, 0.053, 0.025, 0.026]
#'''

# Define weights for each month and weather variable according to Hosseini et al. (2020) for Vancouver
'''
FSWeightsTdryb = [0.786, 0.825, 0.641, 0.595, 0.694, 0.535, 0.429, 0.530, 0.443, 0.718, 0.806, 0.862]
FSWeightsPressure = [0.138, 0.117, 0.243, 0.254, 0.206, 0.394, 0.490, 0.405, 0.434, 0.194, 0.128, 0.098]
FSWeightsRad = [0.046, 0.028, 0.069, 0.073, 0.056, 0.035, 0.039, 0.040, 0.073, 0.058, 0.047, 0.022]
FSWeightsWindS = [0.030, 0.030, 0.047, 0.078, 0.044, 0.036, 0.042, 0.025, 0.050, 0.030, 0.019, 0.018]
'''

# Calculations
########################################################################################################################

# Loading and analyzing EPW Data
########################################################################################################################

# Define matrix for EPW data
EPWPastTdryb = numpy.zeros((NYearsPast, EPWNPoints))
EPWPastPressure = numpy.zeros((NYearsPast, EPWNPoints))
EPWPastRad = numpy.zeros((NYearsPast, EPWNPoints))
EPWPastWindS = numpy.zeros((NYearsPast, EPWNPoints))

# Change the directory
os.chdir(PathEPWHistorical)

# Read text files and load data to matrices
i = 0
for file in os.listdir():
    data = pandas.read_csv(file, header=NSkipHeader)
    data = data.values

    # Note pandas.read_csv() cannot read first hour of data in EPW file
    # Consider the second hour data as first hour data
    EPWPastTdryb[i][0] = CelciusToKelvin + data[0, 6]
    EPWPastPressure[i][0] = data[0, 9]
    # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
    EPWPastRad[i][0] = data[0, 12] + data[0, 14] + data[0, 15]
    EPWPastWindS[i][0] = data[0, 21]

    for l in range(1, EPWNPoints):
        EPWPastTdryb[i][l] = CelciusToKelvin + data[l-1, 6]
        EPWPastPressure[i][l] = data[l-1, 9]
        # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
        EPWPastRad[i][l] = data[l-1, 12] + data[l-1, 14] + data[l-1, 15]
        EPWPastWindS[i][l] = data[l-1, 21]
    i = i + 1

print('\nFinished reading historical EPW files.')
print('Last file: ', file)

# Define matrix for EPW data

EPWValidTdryb = numpy.zeros((NYearsValid, EPWNPoints))
EPWValidPressure = numpy.zeros((NYearsValid, EPWNPoints))
EPWValidRad = numpy.zeros((NYearsValid, EPWNPoints))
EPWValidWindS = numpy.zeros((NYearsValid, EPWNPoints))

# Change the directory
os.chdir(PathEPWValidation)

# Read text files and load data to matrices
i = 0
for file in os.listdir():
    data = pandas.read_csv(file, header=NSkipHeader)
    data = data.values

    # Note pandas.read_csv() cannot read first hour of data in EPW file
    # Consider the second hour data as first hour data
    EPWValidTdryb[i][0] = CelciusToKelvin + data[0, 6]
    EPWValidPressure[i][0] = data[0, 9]
    # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
    EPWValidRad[i][0] = data[0, 12] + data[0, 14] + data[0, 15]
    EPWValidWindS[i][0] = data[0, 21]

    for l in range(1, EPWNPoints):
        EPWValidTdryb[i][l] = CelciusToKelvin + data[l-1, 6]
        EPWValidPressure[i][l] = data[l-1, 9]
        # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
        EPWValidRad[i][l] = data[l-1, 12] + data[l-1, 14] + data[l-1, 15]
        EPWValidWindS[i][l] = data[l-1, 21]
    i = i + 1

print('\nFinished reading validation EPW files.')
print('Last file: ', file)

# Perform block averaging of data to convert 24 hours to daily averages
# Define matrix for EPW daily data
EPWDailyPastTdryb = numpy.zeros((NYearsPast, DailyNPoints))
EPWDailyPastPressure = numpy.zeros((NYearsPast, DailyNPoints))
EPWDailyPastRad = numpy.zeros((NYearsPast, DailyNPoints))
EPWDailyPastWindS = numpy.zeros((NYearsPast, DailyNPoints))

EPWDailyValidTdryb = numpy.zeros((NYearsPast, DailyNPoints))
EPWDailyValidPressure = numpy.zeros((NYearsPast, DailyNPoints))
EPWDailyValidRad = numpy.zeros((NYearsPast, DailyNPoints))
EPWDailyValidWindS = numpy.zeros((NYearsPast, DailyNPoints))

#Past
for i in range(0, NYearsPast):
    for j in range(0, DailyNPoints):
        EPWDailyPastTdryb[i][j] = numpy.nanmean(EPWPastTdryb[i][j*NHoursInDay:j*NHoursInDay+NHoursInDay])
        EPWDailyPastPressure[i][j] = numpy.nanmean(EPWPastPressure[i][j * NHoursInDay:j * NHoursInDay + NHoursInDay])
        EPWDailyPastRad[i][j] = numpy.nanmean(EPWPastRad[i][j * NHoursInDay:j * NHoursInDay + NHoursInDay])
        EPWDailyPastWindS[i][j] = numpy.nanmean(EPWPastWindS[i][j * NHoursInDay:j * NHoursInDay + NHoursInDay])

#Validation
for i in range(0, NYearsValid):
    for j in range(0, DailyNPoints):
        EPWDailyValidTdryb[i][j] = numpy.nanmean(EPWValidTdryb[i][j*NHoursInDay:j*NHoursInDay+NHoursInDay])
        EPWDailyValidPressure[i][j] = numpy.nanmean(EPWValidPressure[i][j*NHoursInDay:j*NHoursInDay+NHoursInDay])
        EPWDailyValidRad[i][j] = numpy.nanmean(EPWValidRad[i][j*NHoursInDay:j*NHoursInDay+NHoursInDay])
        EPWDailyValidWindS[i][j] = numpy.nanmean(EPWValidWindS[i][j*NHoursInDay:j*NHoursInDay+NHoursInDay])

# Cumulative Distribution Function (CDF) Analysis
# Define matrix for quantile data
EPWPastTdrybQuantiles = numpy.zeros((NYearsPast, NMonths, QuantileStepsMax))
EPWPastPressureQuantiles = numpy.zeros((NYearsPast, NMonths, QuantileStepsMax))
EPWPastRadQuantiles = numpy.zeros((NYearsPast, NMonths, QuantileStepsMax))
EPWPastWindSQuantiles = numpy.zeros((NYearsPast, NMonths, QuantileStepsMax))

EPWValidTdrybQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
EPWValidPressureQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
EPWValidRadQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
EPWValidWindSQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))

# Analyze EPW dataset and calculate and store the quantiles
for i in range(0, NYearsPast):
    for j in range(0, NMonths):
        for k in range(0, QuantileSteps[j]):
            EPWPastTdrybQuantiles[i][j][k] = numpy.nanpercentile(EPWDailyPastTdryb[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            EPWPastPressureQuantiles[i][j][k] = numpy.nanpercentile(EPWDailyPastPressure[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            EPWPastRadQuantiles[i][j][k] = numpy.nanpercentile(EPWDailyPastRad[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            EPWPastWindSQuantiles[i][j][k] = numpy.nanpercentile(EPWDailyPastWindS[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])

for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        for k in range(0, QuantileSteps[j]):
            EPWValidTdrybQuantiles[i][j][k] = numpy.nanpercentile(EPWDailyValidTdryb[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            EPWValidPressureQuantiles[i][j][k] = numpy.nanpercentile(EPWDailyValidPressure[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            EPWValidRadQuantiles[i][j][k] = numpy.nanpercentile(EPWDailyValidRad[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            EPWValidWindSQuantiles[i][j][k] = numpy.nanpercentile(EPWDailyValidWindS[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])

print('\nFinished quantile-quantile analysis historical and validation EPW data.')

# Read bias corrected climate data for the validation period

CMValidTdrybQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidPressureQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidRadQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidWindSQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))

# Change the directory
os.chdir(PathCMValidation)

# Read text files and load data to matrices
i = 0
for file in os.listdir():
    # Check whether file is in text format or not
    if file.endswith(".txt"):
        file_path = f"{PathCMValidation}\{file}"
        # call read text file function
        read_text_file(file_path)
    data = numpy.loadtxt(file)

    for l in range(0, DailyNPoints):
        CMValidTdrybQuantiles[i][int(data[l][0])][int(data[l][1])] = data[l][2]
        CMValidPressureQuantiles[i][int(data[l][0])][int(data[l][1])] = data[l][3]
        CMValidRadQuantiles[i][int(data[l][0])][int(data[l][1])] = data[l][4]
        CMValidWindSQuantiles[i][int(data[l][0])][int(data[l][1])] = data[l][5]
    i = i + 1

print('\nFinished reading quantile-quantile validation CM data.')
print('Last file: ', file)

# Calculate Finkelstein-Schafer (FS) Method for validation period
# Define matrix for FS for each validation year and month
FSValidTdryb = numpy.zeros((NYearsValid, NMonths, NYearsPast, NMonths))
FSValidPressure = numpy.zeros((NYearsValid, NMonths, NYearsPast, NMonths))
FSValidRad = numpy.zeros((NYearsValid, NMonths, NYearsPast, NMonths))
FSValidWindS = numpy.zeros((NYearsValid, NMonths, NYearsPast, NMonths))
WSValid = numpy.zeros((NYearsValid, NMonths, NYearsPast, NMonths))

for i in range(0, NYearsValid):
        for j in range(0, NMonths):
            for k in range(0, NYearsPast):
                for l in range(0, NMonths):
                    FSValidTdryb[i][j][k][l] = numpy.nansum(
                        numpy.abs(CMValidTdrybQuantiles[i][j][:] - EPWPastTdrybQuantiles[k][l][:])) / QuantileSteps[j]
                    FSValidPressure[i][j][k][l] = numpy.nansum(
                        numpy.abs(CMValidPressureQuantiles[i][j][:] - EPWPastPressureQuantiles[k][l][:])) / QuantileSteps[j]
                    FSValidRad[i][j][k][l] = numpy.nansum(
                        numpy.abs(CMValidRadQuantiles[i][j][:] - EPWPastRadQuantiles[k][l][:])) / QuantileSteps[j]
                    FSValidWindS[i][j][k][l] = numpy.nansum(
                        numpy.abs(CMValidWindSQuantiles[i][j][:] - EPWPastWindSQuantiles[k][l][:])) / QuantileSteps[j]
                    WSValid[i][j][k][l] = FSWeightsTdryb[j] * FSValidTdryb[i][j][k][l] + FSWeightsPressure[j] * FSValidPressure[i][j][k][l] + \
                       FSWeightsRad[j] * FSValidRad[i][j][k][l] + FSWeightsWindS[j] * FSValidWindS[i][j][k][l]

print('\nApplied Finkelstein-Schafer (FS) method for validation CM data.')

# Look for minimum weighted sum for every month and year of the validation period
# Check within all the past months to find the lowest weighted sum
ValidIndex = numpy.zeros((NYearsValid, NMonths, 2))

for i in range(0, NYearsValid):
        for j in range(0, NMonths):
            MinimumWS = 10e6
            for k in range(0, NYearsPast):
                for l in range(0, NMonths):
                    if WSValid[i][j][k][l] < MinimumWS:
                        MinimumWS = WSValid[i][j][k][l]
                        ValidIndex[i][j][0] = k
                        ValidIndex[i][j][1] = l

print('\nCalculated year/month indices for matching EPW records for validation CM data.')

print('\nSample: best matching year index: ', int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][0]))
print('Sample: best matching month index: ', int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][1]))

CDFDaily = numpy.zeros((NMonths, QuantileStepsMax))
for i in range(0, NMonths):
    for j in range(1, QuantileSteps[i] + 1):
        CDFDaily[i][j-1] = j/QuantileSteps[i]

print('\nPlot sample results ...')

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
for i in range(0, NYearsPast):
    for j in range(0, NMonths):
        plt.plot(EPWPastTdrybQuantiles[i][j][0:QuantileSteps[j]]-CelciusToKelvin, CDFDaily[j][0:QuantileSteps[j]], 'k-', alpha = 0.1)
plt.plot(CMValidTdrybQuantiles[PlotYearValidIndex][PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]]-CelciusToKelvin, CDFDaily[PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], 'b-', linewidth = 3)
plt.plot(EPWPastTdrybQuantiles[int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][0])][int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][1])][0:QuantileSteps[PlotMonthValidIndex]]-CelciusToKelvin, CDFDaily[PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], 'r-', linewidth = 3)
plt.xlabel('Dry Bulb Temperature [$^{\circ}$C]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
EPWPastAllTdryb = mpatches.Patch(color='k', label='All Possibilities', alpha = 0.5)
CMofInterestTdryb = mpatches.Patch(color='blue', label='Interest')
EPWBestMatchTdryb = mpatches.Patch(color='red', label='Best Match')
plt.legend(handles=[EPWPastAllTdryb, CMofInterestTdryb, EPWBestMatchTdryb], fontsize=20)
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQTdryb.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
for i in range(0, NYearsPast):
    for j in range(0, NMonths):
        plt.plot(EPWPastPressureQuantiles[i][j][0:QuantileSteps[j]], CDFDaily[j][0:QuantileSteps[j]], 'k-', alpha = 0.1)
plt.plot(CMValidPressureQuantiles[PlotYearValidIndex][PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], CDFDaily[PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], 'b-', linewidth = 3)
plt.plot(EPWPastPressureQuantiles[int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][0])][int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][1])][0:QuantileSteps[PlotMonthValidIndex]], CDFDaily[PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], 'r-', linewidth = 3)
plt.xlabel('Pressure [Pa]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
EPWPastAllPs = mpatches.Patch(color='k', label='All Possibilities', alpha = 0.5)
CMofInterestPs = mpatches.Patch(color='blue', label='Interest')
EPWBestMatchPs = mpatches.Patch(color='red', label='Best Match')
plt.legend(handles=[EPWPastAllPs, CMofInterestPs, EPWBestMatchPs], fontsize=20)
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQPressure.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
for i in range(0, NYearsPast):
    for j in range(0, NMonths):
        plt.plot(EPWPastRadQuantiles[i][j][0:QuantileSteps[j]], CDFDaily[j][0:QuantileSteps[j]], 'k-', alpha = 0.1)
plt.plot(CMValidRadQuantiles[PlotYearValidIndex][PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], CDFDaily[PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], 'b-', linewidth = 3)
plt.plot(EPWPastRadQuantiles[int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][0])][int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][1])][0:QuantileSteps[PlotMonthValidIndex]], CDFDaily[PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], 'r-', linewidth = 3)
plt.xlabel('Global Horiz. Rad. [W m$^{-2}$]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
EPWPastAllRad = mpatches.Patch(color='k', label='All Possibilities', alpha = 0.5)
CMofInterestRad = mpatches.Patch(color='blue', label='Interest')
EPWBestMatchRad = mpatches.Patch(color='red', label='Best Match')
plt.legend(handles=[EPWPastAllRad, CMofInterestRad, EPWBestMatchRad], fontsize=20)
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQRad.pdf', dpi=600)
fig.show()

fig = plt.figure(figsize=(8,6))
for i in range(0, NYearsPast):
    for j in range(0, NMonths):
        plt.plot(EPWPastWindSQuantiles[i][j][0:QuantileSteps[j]], CDFDaily[j][0:QuantileSteps[j]], 'k-', alpha = 0.1)
plt.plot(CMValidWindSQuantiles[PlotYearValidIndex][PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], CDFDaily[PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], 'b-', linewidth = 3)
plt.plot(EPWPastWindSQuantiles[int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][0])][int(ValidIndex[PlotYearValidIndex][PlotMonthValidIndex][1])][0:QuantileSteps[PlotMonthValidIndex]], CDFDaily[PlotMonthValidIndex][0:QuantileSteps[PlotMonthValidIndex]], 'r-', linewidth = 3)
plt.xlabel('Wind Speed [m s$^{-1}$]', fontsize=20)
plt.ylabel('CDF [-]', fontsize=20)
EPWPastAllWindS = mpatches.Patch(color='k', label='All Possibilities', alpha = 0.5)
CMofInterestWindS = mpatches.Patch(color='blue', label='Interest')
EPWBestMatchWindS = mpatches.Patch(color='red', label='Best Match')
plt.legend(handles=[EPWPastAllWindS, CMofInterestWindS, EPWBestMatchWindS], fontsize=20)
plt.ylim((0,1))
plt.tight_layout()
plt.savefig('QQWindS.pdf', dpi=600)
fig.show()

plt.show()