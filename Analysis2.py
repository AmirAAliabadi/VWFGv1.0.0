# Quantile-Quantile analysis of EPW and CM to correct for bias
# Output bias of non-corrected CM vs EPW data for the historical period
# Output bias of non-corrected CM vs EPW data for the validation period
# Output bias of corrected CM vs EPW data for the validation period
# Output daily-averaged mean/std of EPW data for the past period
# Output daily-averaged mean/std of non-corrected CM data for the validation period
# Output daily-averaged mean/std of corrected CM data for the validation period
# Atmospheric Innovations Research (AIR) Laboratory, University of Guelph
# Developed by Amir A. Aliabadi and Rachel M. McLeod
# Last updated: 2022-03-16

import numpy
import pandas
import os
import matplotlib.pyplot as plt

# Function definitions
########################################################################################################################

def read_text_file(file_path):
    with open(file_path, 'r') as f:
        print(f.read())

# Constants of simulation
########################################################################################################################

# Historical period: 1980-1999
# Validation period: 2007-2020
# Future period: 2021-2100

# Define directories and file names
PathEPWHistorical = ("C:/GoogleDrive/U Guelph/Projects/VWFG/EPW/EPW-Historical")
PathEPWValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/EPW/EPW-Validation")
PathCMHistorical = ("C:/GoogleDrive/U Guelph/Projects/VWFG/CM/CM-Historical")
PathCMValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/CM/CM-Validation")
PathCMFuture = ("C:/GoogleDrive/U Guelph/Projects/VWFG/CM/CM-Future")
PathResults = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/")
PathResultsCMValidationBiasCorrected = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/CM-Validation-Bias-Corrected")
PathResultsCMFutureBiasCorrected = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/CM-Future-Bias-Corrected")

FileNameBiasCMEPWPast = PathResults + 'Bias-CM-EPW-Past.txt'
FileNameBiasCMEPWValidNoCorrection = PathResults + 'Bias-CM-EPW-Validation-No-Correction.txt'
FileNameBiasCMEPWValidCorrection = PathResults + 'Bias-CM-EPW-Validation-Correction.txt'
FileNameMeanStdEPWPast = PathResults + 'Mean-Std-EPW-Past.txt'
FileNameMeanStdCMValidCorrection = PathResults + 'Mean-Std-CM-Validation-Correction.txt'
FileNameMeanStdCMFutureCorrection = PathResults + 'Mean-Std-CM-Future-Correction.txt'

# Specify the number of past historical, validation, and future years in the dataset
# Specify the number of months in each year
NYearsPast = 20
NYearsValid = 14
NYearsFuture = 80
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
    EPWValidWindS[0][0] = data[0, 21]

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

# Define matrices for quantile data combining multiple years for the past
EPWPastMultiYearTdrybQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
EPWPastMultiYearPressureQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
EPWPastMultiYearRadQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
EPWPastMultiYearWindSQuantiles = numpy.zeros((NMonths, QuantileStepsMax))

# Calculate CFD of the past for multi year period
for j in range(0, NMonths):
    for k in range(0, QuantileSteps[j]):
        SumTdryb = 0
        SumPressure = 0
        SumRad = 0
        SumWindS = 0
        for i in range(0, NYearsPast):
            SumTdryb = SumTdryb + EPWPastTdrybQuantiles[i][j][k]
            SumPressure = SumPressure + EPWPastPressureQuantiles[i][j][k]
            SumRad = SumRad + EPWPastRadQuantiles[i][j][k]
            SumWindS = SumWindS + EPWPastWindSQuantiles[i][j][k]

        EPWPastMultiYearTdrybQuantiles[j][k] = SumTdryb / NYearsPast
        EPWPastMultiYearPressureQuantiles[j][k] = SumPressure / NYearsPast
        EPWPastMultiYearRadQuantiles[j][k] = SumRad / NYearsPast
        EPWPastMultiYearWindSQuantiles[j][k] = SumWindS / NYearsPast

# Define matrices for quantile data combining multiple years for validation period
EPWValidMultiYearTdrybQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
EPWValidMultiYearPressureQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
EPWValidMultiYearRadQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
EPWValidMultiYearWindSQuantiles = numpy.zeros((NMonths, QuantileStepsMax))

# Calculate CFD of the validation period for multi year period
for j in range(0, NMonths):
    for k in range(0, QuantileSteps[j]):
        SumTdryb = 0
        SumPressure = 0
        SumRad = 0
        SumWindS = 0
        for i in range(0, NYearsValid):
            SumTdryb = SumTdryb + EPWValidTdrybQuantiles[i][j][k]
            SumPressure = SumPressure + EPWValidPressureQuantiles[i][j][k]
            SumRad = SumRad + EPWValidRadQuantiles[i][j][k]
            SumWindS = SumWindS + EPWValidWindSQuantiles[i][j][k]

        EPWValidMultiYearTdrybQuantiles[j][k] = SumTdryb / NYearsValid
        EPWValidMultiYearPressureQuantiles[j][k] = SumPressure / NYearsValid
        EPWValidMultiYearRadQuantiles[j][k] = SumRad / NYearsValid
        EPWValidMultiYearWindSQuantiles[j][k] = SumWindS / NYearsValid

# Loading and Analyzing Climate Model Data
########################################################################################################################

# Read climate data from climate model for the past

CMDailyPastTdryb = numpy.zeros((NYearsPast, DailyNPoints))
CMDailyPastPressure = numpy.zeros((NYearsPast, DailyNPoints))
CMDailyPastRad = numpy.zeros((NYearsPast, DailyNPoints))
CMDailyPastWindS = numpy.zeros((NYearsPast, DailyNPoints))

# Change the directory
os.chdir(PathCMHistorical)

# Read text files and load data to matrices
i = 0
for file in os.listdir():
    # Check whether file is in text format or not
    if file.endswith(".txt"):
        file_path = f"{PathCMHistorical}\{file}"
        # call read text file function
        read_text_file(file_path)
    data = numpy.loadtxt(file)
    CMDailyPastTdryb[i][:] = data[:, 3]
    CMDailyPastPressure[i][:] = data[:, 5]
    CMDailyPastRad[i][:] = data[:, 6]
    CMDailyPastWindS[i][:] = data[:, 4]
    i = i + 1

print('\nFinished reading historical CM files.')
print('Last file: ', file)

# Read climate data from climate model for validation

CMDailyValidTdryb = numpy.zeros((NYearsValid, DailyNPoints))
CMDailyValidPressure = numpy.zeros((NYearsValid, DailyNPoints))
CMDailyValidRad = numpy.zeros((NYearsValid, DailyNPoints))
CMDailyValidWindS = numpy.zeros((NYearsValid, DailyNPoints))

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
    CMDailyValidTdryb[i][:] = data[:, 3]
    CMDailyValidPressure[i][:] = data[:, 5]
    CMDailyValidRad[i][:] = data[:, 6]
    CMDailyValidWindS[i][:] = data[:, 4]
    i = i + 1

print('\nFinished reading validation CM files.')
print('Last file: ', file)

# Read climate data from climate model for future

CMDailyFutureTdryb = numpy.zeros((NYearsFuture, DailyNPoints))
CMDailyFuturePressure = numpy.zeros((NYearsFuture, DailyNPoints))
CMDailyFutureRad = numpy.zeros((NYearsFuture, DailyNPoints))
CMDailyFutureWindS = numpy.zeros((NYearsFuture, DailyNPoints))

# Change the directory
os.chdir(PathCMFuture)

# Read text files and load data to matrices
i = 0
for file in os.listdir():
    # Check whether file is in text format or not
    if file.endswith(".txt"):
        file_path = f"{PathCMFuture}\{file}"
        # call read text file function
        read_text_file(file_path)
    data = numpy.loadtxt(file)
    CMDailyFutureTdryb[i][:] = data[:, 3]
    CMDailyFuturePressure[i][:] = data[:, 5]
    CMDailyFutureRad[i][:] = data[:, 6]
    CMDailyFutureWindS[i][:] = data[:, 4]
    i = i + 1

print('\nFinished reading future CM files.')
print('Last file: ', file)

# Cumulative Distribution Function (CDF) Analysis
# Define matrix for quantile data
CMPastTdrybQuantiles = numpy.zeros((NYearsPast, NMonths, QuantileStepsMax))
CMPastPressureQuantiles = numpy.zeros((NYearsPast, NMonths, QuantileStepsMax))
CMPastRadQuantiles = numpy.zeros((NYearsPast, NMonths, QuantileStepsMax))
CMPastWindSQuantiles = numpy.zeros((NYearsPast, NMonths, QuantileStepsMax))

CMValidTdrybQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidPressureQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidRadQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidWindSQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))

CMFutureTdrybQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFuturePressureQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFutureRadQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFutureWindSQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))

# Analyze CM dataset and calculate and store the quantiles
for i in range(0, NYearsPast):
    for j in range(0, NMonths):
        for k in range(0, QuantileSteps[j]):
            CMPastTdrybQuantiles[i][j][k] = numpy.nanpercentile(CMDailyPastTdryb[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMPastPressureQuantiles[i][j][k] = numpy.nanpercentile(CMDailyPastPressure[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMPastRadQuantiles[i][j][k] = numpy.nanpercentile(CMDailyPastRad[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMPastWindSQuantiles[i][j][k] = numpy.nanpercentile(CMDailyPastWindS[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])

for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        for k in range(0, QuantileSteps[j]):
            CMValidTdrybQuantiles[i][j][k] = numpy.nanpercentile(CMDailyValidTdryb[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMValidPressureQuantiles[i][j][k] = numpy.nanpercentile(CMDailyValidPressure[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMValidRadQuantiles[i][j][k] = numpy.nanpercentile(CMDailyValidRad[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMValidWindSQuantiles[i][j][k] = numpy.nanpercentile(CMDailyValidWindS[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])

for i in range(0, NYearsFuture):
    for j in range(0, NMonths):
        for k in range(0, QuantileSteps[j]):
            CMFutureTdrybQuantiles[i][j][k] = numpy.nanpercentile(CMDailyFutureTdryb[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMFuturePressureQuantiles[i][j][k] = numpy.nanpercentile(CMDailyFuturePressure[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMFutureRadQuantiles[i][j][k] = numpy.nanpercentile(CMDailyFutureRad[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])
            CMFutureWindSQuantiles[i][j][k] = numpy.nanpercentile(CMDailyFutureWindS[i][MonthIndex[j]:MonthIndex[j]+QuantileSteps[j]],
                                                             100 * (k + 1) / QuantileSteps[j])

print('\nFinished quantile-quantile analysis historical, validation, and future CM data.')

# Define matrices for quantile data combining multiple years for the past

CMPastMultiYearTdrybQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMPastMultiYearPressureQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMPastMultiYearRadQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMPastMultiYearWindSQuantiles = numpy.zeros((NMonths, QuantileStepsMax))

# Calculate CFD of the past for multi year period
for j in range(0, NMonths):
    for k in range(0, QuantileSteps[j]):
        SumTdryb = 0
        SumPressure = 0
        SumRad = 0
        SumWindS = 0
        for i in range(0, NYearsPast):
            SumTdryb = SumTdryb + CMPastTdrybQuantiles[i][j][k]
            SumPressure = SumPressure + CMPastPressureQuantiles[i][j][k]
            SumRad = SumRad + CMPastRadQuantiles[i][j][k]
            SumWindS = SumWindS + CMPastWindSQuantiles[i][j][k]

        CMPastMultiYearTdrybQuantiles[j][k] = SumTdryb / NYearsPast
        CMPastMultiYearPressureQuantiles[j][k] = SumPressure / NYearsPast
        CMPastMultiYearRadQuantiles[j][k] = SumRad / NYearsPast
        CMPastMultiYearWindSQuantiles[j][k] = SumWindS / NYearsPast

# Define matrices for quantile data combining multiple years for the validation period

CMValidMultiYearTdrybQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMValidMultiYearPressureQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMValidMultiYearRadQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMValidMultiYearWindSQuantiles = numpy.zeros((NMonths, QuantileStepsMax))

# Calculate CFD of the validation period for multi year period
for j in range(0, NMonths):
    for k in range(0, QuantileSteps[j]):
        SumTdryb = 0
        SumPressure = 0
        SumRad = 0
        SumWindS = 0
        for i in range(0, NYearsValid):
            SumTdryb = SumTdryb + CMValidTdrybQuantiles[i][j][k]
            SumPressure = SumPressure + CMValidPressureQuantiles[i][j][k]
            SumRad = SumRad + CMValidRadQuantiles[i][j][k]
            SumWindS = SumWindS + CMValidWindSQuantiles[i][j][k]

        CMValidMultiYearTdrybQuantiles[j][k] = SumTdryb / NYearsValid
        CMValidMultiYearPressureQuantiles[j][k] = SumPressure / NYearsValid
        CMValidMultiYearRadQuantiles[j][k] = SumRad / NYearsValid
        CMValidMultiYearWindSQuantiles[j][k] = SumWindS / NYearsValid

print('\nFinished multi-year quantile-quantile analysis of historical and validation CM data.')

# Define matrices for biases
BiasPastTdryb = numpy.zeros((NMonths, QuantileStepsMax))
BiasPastPressure = numpy.zeros((NMonths, QuantileStepsMax))
BiasPastRad = numpy.zeros((NMonths, QuantileStepsMax))
BiasPastWindS = numpy.zeros((NMonths, QuantileStepsMax))

BiasValidNoCorrectionTdryb = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidNoCorrectionPressure = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidNoCorrectionRad = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidNoCorrectionWindS = numpy.zeros((NMonths, QuantileStepsMax))

BiasValidCorrectionTdryb = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidCorrectionPressure = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidCorrectionRad = numpy.zeros((NMonths, QuantileStepsMax))
BiasValidCorrectionWindS = numpy.zeros((NMonths, QuantileStepsMax))

for j in range(0, NMonths):
    for k in range(0, QuantileSteps[j]):
        BiasPastTdryb[j][k] = CMPastMultiYearTdrybQuantiles[j][k] - EPWPastMultiYearTdrybQuantiles[j][k]
        BiasPastPressure[j][k] = CMPastMultiYearPressureQuantiles[j][k] - EPWPastMultiYearPressureQuantiles[j][k]
        BiasPastRad[j][k] = CMPastMultiYearRadQuantiles[j][k] - EPWPastMultiYearRadQuantiles[j][k]
        BiasPastWindS[j][k] = CMPastMultiYearWindSQuantiles[j][k] - EPWPastMultiYearWindSQuantiles[j][k]

        BiasValidNoCorrectionTdryb[j][k] = CMValidMultiYearTdrybQuantiles[j][k] - EPWValidMultiYearTdrybQuantiles[j][k]
        BiasValidNoCorrectionPressure[j][k] = CMValidMultiYearPressureQuantiles[j][k] - EPWValidMultiYearPressureQuantiles[j][k]
        BiasValidNoCorrectionRad[j][k] = CMValidMultiYearRadQuantiles[j][k] - EPWValidMultiYearRadQuantiles[j][k]
        BiasValidNoCorrectionWindS[j][k] = CMValidMultiYearWindSQuantiles[j][k] - EPWValidMultiYearWindSQuantiles[j][k]

# Subtract bias from CM data for the validation and future period and write to new files

CMValidCorrectedTdrybQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidCorrectedPressureQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidCorrectedRadQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))
CMValidCorrectedWindSQuantiles = numpy.zeros((NYearsValid, NMonths, QuantileStepsMax))

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

    outputFileName = file + '-CDF-' + 'Bias-Corrected'

    outputFile = open(PathResultsCMValidationBiasCorrected+'/'+outputFileName, "w")
    outputFile.write(
        "# 0: Month [-] \t 1: Quantile [-] \t 2: Tdryb [K] \t 3: Pressure [Pa] \t 4: Rad [W m-2] \t 5: WindS [m s-1] \n")

    for j in range(0, NMonths):
        for k in range(0, QuantileSteps[j]):
            CMValidCorrectedTdrybQuantiles[i][j][k] = CMValidTdrybQuantiles[i][j][k] - BiasPastTdryb[j][k]
            CMValidCorrectedPressureQuantiles[i][j][k] = CMValidPressureQuantiles[i][j][k] - BiasPastPressure[j][k]
            CMValidCorrectedRadQuantiles[i][j][k] = CMValidRadQuantiles[i][j][k] - BiasPastRad[j][k]
            CMValidCorrectedWindSQuantiles[i][j][k] = CMValidWindSQuantiles[i][j][k] - BiasPastWindS[j][k]

            outputFile.write("%i \t %i \t %f \t %f \t %f \t %f \n"
                                    % (j, k, CMValidCorrectedTdrybQuantiles[i][j][k],
                                       CMValidCorrectedPressureQuantiles[i][j][k],
                                       CMValidCorrectedRadQuantiles[i][j][k],
                                       CMValidCorrectedWindSQuantiles[i][j][k]))

    outputFile.close()
    i = i + 1

print('\nFinished bias correction of validation CM data; Wrote data to file.')
print('Last file: ', outputFileName)

CMFutureCorrectedTdrybQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFutureCorrectedPressureQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFutureCorrectedRadQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFutureCorrectedWindSQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))

# Change the directory
os.chdir(PathCMFuture)

# Read text files and load data to matrices
i = 0
for file in os.listdir():
    # Check whether file is in text format or not
    if file.endswith(".txt"):
        file_path = f"{PathCMFuture}\{file}"
        # call read text file function
        read_text_file(file_path)

    outputFileName = file + '-CDF-' + 'Bias-Corrected'

    outputFile = open(PathResultsCMFutureBiasCorrected+'/'+outputFileName, "w")
    outputFile.write(
        "# 0: Month [-] \t 1: Quantile [-] \t 2: Tdryb [K] \t 3: Pressure [Pa] \t 4: Rad [W m-2] \t 5: WindS [m s-1] \n")

    for j in range(0, NMonths):
        for k in range(0, QuantileSteps[j]):
            CMFutureCorrectedTdrybQuantiles[i][j][k] = CMFutureTdrybQuantiles[i][j][k] - BiasPastTdryb[j][k]
            CMFutureCorrectedPressureQuantiles[i][j][k] = CMFuturePressureQuantiles[i][j][k] - BiasPastPressure[j][k]
            CMFutureCorrectedRadQuantiles[i][j][k] = CMFutureRadQuantiles[i][j][k] - BiasPastRad[j][k]
            CMFutureCorrectedWindSQuantiles[i][j][k] = CMFutureWindSQuantiles[i][j][k] - BiasPastWindS[j][k]

            outputFile.write("%i \t %i \t %f \t %f \t %f \t %f \n"
                                    % (j, k, CMFutureCorrectedTdrybQuantiles[i][j][k],
                                       CMFutureCorrectedPressureQuantiles[i][j][k],
                                       CMFutureCorrectedRadQuantiles[i][j][k],
                                       CMFutureCorrectedWindSQuantiles[i][j][k]))

    outputFile.close()

    i = i + 1

print('\nFinished bias correction of future CM data; Wrote data to file.')
print('Last file: ', outputFileName)

# Define matrices for quantile data combining multiple years for the validation period using corrected data

CMValidCorrectedMultiYearTdrybQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMValidCorrectedMultiYearPressureQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMValidCorrectedMultiYearRadQuantiles = numpy.zeros((NMonths, QuantileStepsMax))
CMValidCorrectedMultiYearWindSQuantiles = numpy.zeros((NMonths, QuantileStepsMax))

# Calculate CFD of the validation period for multi year period
for j in range(0, NMonths):
    for k in range(0, QuantileSteps[j]):
        SumTdryb = 0
        SumPressure = 0
        SumRad = 0
        SumWindS = 0
        for i in range(0, NYearsValid):
            SumTdryb = SumTdryb + CMValidCorrectedTdrybQuantiles[i][j][k]
            SumPressure = SumPressure + CMValidCorrectedPressureQuantiles[i][j][k]
            SumRad = SumRad + CMValidCorrectedRadQuantiles[i][j][k]
            SumWindS = SumWindS + CMValidCorrectedWindSQuantiles[i][j][k]

        CMValidCorrectedMultiYearTdrybQuantiles[j][k] = SumTdryb / NYearsValid
        CMValidCorrectedMultiYearPressureQuantiles[j][k] = SumPressure / NYearsValid
        CMValidCorrectedMultiYearRadQuantiles[j][k] = SumRad / NYearsValid
        CMValidCorrectedMultiYearWindSQuantiles[j][k] = SumWindS / NYearsValid

# Calculate bias after correction
for j in range(0, NMonths):            
    for k in range(0, QuantileSteps[j]):
        BiasValidCorrectionTdryb[j][k] = CMValidCorrectedMultiYearTdrybQuantiles[j][k] - EPWValidMultiYearTdrybQuantiles[j][k]
        BiasValidCorrectionPressure[j][k] = CMValidCorrectedMultiYearPressureQuantiles[j][k] - EPWValidMultiYearPressureQuantiles[j][k]
        BiasValidCorrectionRad[j][k] = CMValidCorrectedMultiYearRadQuantiles[j][k] - EPWValidMultiYearRadQuantiles[j][k]
        BiasValidCorrectionWindS[j][k] = CMValidCorrectedMultiYearWindSQuantiles[j][k] - EPWValidMultiYearWindSQuantiles[j][k]

# Writing results to file
FileBiasCMEPWPast = open(FileNameBiasCMEPWPast, "w")
FileBiasCMEPWPast.write("# 0: Month [-] \t 1: Quantile [-] \t 2: Tdryb [K] \t 3: Pressure [Pa] \t 4: Rad [W m-2] \t 5: WindS [m s-1] \n")

FileBiasCMEPWValidNoCorrection = open(FileNameBiasCMEPWValidNoCorrection, "w")
FileBiasCMEPWValidNoCorrection.write("# 0: Month [-] \t 1: Quantile [-] \t 2: Tdryb [K] \t 3: Pressure [Pa] \t 4: Rad [W m-2] \t 5: WindS [m s-1] \n")

FileBiasCMEPWValidCorrection = open(FileNameBiasCMEPWValidCorrection, "w")
FileBiasCMEPWValidCorrection.write("# 0: Month [-] \t 1: Quantile [-] \t 2: Tdryb [K] \t 3: Pressure [Pa] \t 4: Rad [W m-2] \t 5: WindS [m s-1] \n")

for j in range(0, NMonths):
    for k in range(0, QuantileSteps[j]):
        FileBiasCMEPWPast.write("%i \t %i \t %f \t %f \t %f \t %f \n"
                          % (j, k, BiasPastTdryb[j][k], BiasPastPressure[j][k], BiasPastRad[j][k], BiasPastWindS[j][k]))
        FileBiasCMEPWValidNoCorrection.write("%i \t %i \t %f \t %f \t %f \t %f \n"
                          % (j, k, BiasValidNoCorrectionTdryb[j][k], BiasValidNoCorrectionPressure[j][k], BiasValidNoCorrectionRad[j][k],
                                   BiasValidNoCorrectionWindS[j][k]))
        FileBiasCMEPWValidCorrection.write("%i \t %i \t %f \t %f \t %f \t %f \n"
                          % (j, k, BiasValidCorrectionTdryb[j][k], BiasValidCorrectionPressure[j][k], BiasValidCorrectionRad[j][k],
                                   BiasValidCorrectionWindS[j][k]))

FileBiasCMEPWPast.close()
FileBiasCMEPWValidNoCorrection.close()
FileBiasCMEPWValidCorrection.close()

print('\nFinished bias calculation of historical CM data.')
print('Wrote data to file: ', FileNameBiasCMEPWPast)

print('\nFinished bias calculation of validation CM data without bias correction.')
print('Wrote data to file: ', FileNameBiasCMEPWValidNoCorrection)

print('\nFinished bias calculation of validation CM data with bias correction.')
print('Wrote data to file: ', FileNameBiasCMEPWValidCorrection)

# Calculate mean and standard deviation of daily weather variables from bias-corrected climate model for each month of the validation period
CMValidCorrectedMeanTdryb = numpy.zeros((NYearsValid, NMonths))
CMValidCorrectedMeanPressure = numpy.zeros((NYearsValid, NMonths))
CMValidCorrectedMeanRad = numpy.zeros((NYearsValid, NMonths))
CMValidCorrectedMeanWindS = numpy.zeros((NYearsValid, NMonths))

CMValidCorrectedStdTdryb = numpy.zeros((NYearsValid, NMonths))
CMValidCorrectedStdPressure = numpy.zeros((NYearsValid, NMonths))
CMValidCorrectedStdRad = numpy.zeros((NYearsValid, NMonths))
CMValidCorrectedStdWindS = numpy.zeros((NYearsValid, NMonths))

for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        CMValidCorrectedMeanTdryb[i][j] = numpy.sum(CMValidCorrectedTdrybQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        CMValidCorrectedMeanPressure[i][j] = numpy.sum(CMValidCorrectedPressureQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        CMValidCorrectedMeanRad[i][j] = numpy.sum(CMValidCorrectedRadQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        CMValidCorrectedMeanWindS[i][j] = numpy.sum(CMValidCorrectedWindSQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]

        CMValidCorrectedStdTdryb[i][j] = numpy.sqrt(numpy.sum(
            (CMValidCorrectedTdrybQuantiles[i][j][0:QuantileSteps[j]] - CMValidCorrectedMeanTdryb[i][j]) ** 2) /
                                                    QuantileSteps[j])
        CMValidCorrectedStdPressure[i][j] = numpy.sqrt(numpy.sum(
            (CMValidCorrectedPressureQuantiles[i][j][0:QuantileSteps[j]] - CMValidCorrectedMeanPressure[i][j]) ** 2) /
                                                    QuantileSteps[j])
        CMValidCorrectedStdRad[i][j] = numpy.sqrt(numpy.sum(
            (CMValidCorrectedRadQuantiles[i][j][0:QuantileSteps[j]] - CMValidCorrectedMeanRad[i][j]) ** 2) /
                                                    QuantileSteps[j])
        CMValidCorrectedStdWindS[i][j] = numpy.sqrt(numpy.sum(
            (CMValidCorrectedWindSQuantiles[i][j][0:QuantileSteps[j]] - CMValidCorrectedMeanWindS[i][j]) ** 2) /
                                                    QuantileSteps[j])

# Save data to file
FileMeanStdCMValidCorrection = open(FileNameMeanStdCMValidCorrection, "w")
FileMeanStdCMValidCorrection.write("# 0: Year [-] \t 1: Month [-] \t 2: MeanTdryb [K] \t 3: MeanPressure [Pa] \t 4: MeanRad [W m-2] \t 5: MeanWindS [m s-1]"
                                   "\t 6: StdTdryb [K] \t 7: StdPressure [Pa] \t 8: StdRad [W m-2] \t 9: StdWindS [m s-1] \n")

for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        FileMeanStdCMValidCorrection.write("%i \t %i \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n"
                                % (i, j, CMValidCorrectedMeanTdryb[i][j], CMValidCorrectedMeanPressure[i][j],
                                   CMValidCorrectedMeanRad[i][j], CMValidCorrectedMeanWindS[i][j],
                                   CMValidCorrectedStdTdryb[i][j], CMValidCorrectedStdPressure[i][j],
                                   CMValidCorrectedStdRad[i][j], CMValidCorrectedStdWindS[i][j]))

FileMeanStdCMValidCorrection.close()

print('\nFinished calculation of monthly mean/std of variables for bias corrected validation CM data.')
print('Wrote data to file: ', FileNameMeanStdCMValidCorrection)

CMFutureCorrectedMeanTdryb = numpy.zeros((NYearsFuture, NMonths))
CMFutureCorrectedMeanPressure = numpy.zeros((NYearsFuture, NMonths))
CMFutureCorrectedMeanRad = numpy.zeros((NYearsFuture, NMonths))
CMFutureCorrectedMeanWindS = numpy.zeros((NYearsFuture, NMonths))

CMFutureCorrectedStdTdryb = numpy.zeros((NYearsFuture, NMonths))
CMFutureCorrectedStdPressure = numpy.zeros((NYearsFuture, NMonths))
CMFutureCorrectedStdRad = numpy.zeros((NYearsFuture, NMonths))
CMFutureCorrectedStdWindS = numpy.zeros((NYearsFuture, NMonths))

for i in range(0, NYearsFuture):
    for j in range(0, NMonths):
        CMFutureCorrectedMeanTdryb[i][j] = numpy.sum(CMFutureCorrectedTdrybQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        CMFutureCorrectedMeanPressure[i][j] = numpy.sum(CMFutureCorrectedPressureQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        CMFutureCorrectedMeanRad[i][j] = numpy.sum(CMFutureCorrectedRadQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        CMFutureCorrectedMeanWindS[i][j] = numpy.sum(CMFutureCorrectedWindSQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]

        CMFutureCorrectedStdTdryb[i][j] = numpy.sqrt(numpy.sum(
            (CMFutureCorrectedTdrybQuantiles[i][j][0:QuantileSteps[j]] - CMFutureCorrectedMeanTdryb[i][j]) ** 2) /
                                                    QuantileSteps[j])
        CMFutureCorrectedStdPressure[i][j] = numpy.sqrt(numpy.sum(
            (CMFutureCorrectedPressureQuantiles[i][j][0:QuantileSteps[j]] - CMFutureCorrectedMeanPressure[i][j]) ** 2) /
                                                    QuantileSteps[j])
        CMFutureCorrectedStdRad[i][j] = numpy.sqrt(numpy.sum(
            (CMFutureCorrectedRadQuantiles[i][j][0:QuantileSteps[j]] - CMFutureCorrectedMeanRad[i][j]) ** 2) /
                                                    QuantileSteps[j])
        CMFutureCorrectedStdWindS[i][j] = numpy.sqrt(numpy.sum(
            (CMFutureCorrectedWindSQuantiles[i][j][0:QuantileSteps[j]] - CMFutureCorrectedMeanWindS[i][j]) ** 2) /
                                                    QuantileSteps[j])

# Save data to file
FileMeanStdCMFutureCorrection = open(FileNameMeanStdCMFutureCorrection, "w")
FileMeanStdCMFutureCorrection.write("# 0: Year [-] \t 1: Month [-] \t 2: MeanTdryb [K] \t 3: MeanPressure [Pa] \t 4: MeanRad [W m-2] \t 5: MeanWindS [m s-1]"
                                   "\t 6: StdTdryb [K] \t 7: StdPressure [Pa] \t 8: StdRad [W m-2] \t 9: StdWindS [m s-1] \n")

for i in range(0, NYearsFuture):
    for j in range(0, NMonths):
        FileMeanStdCMFutureCorrection.write("%i \t %i \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n"
                                % (i, j, CMFutureCorrectedMeanTdryb[i][j], CMFutureCorrectedMeanPressure[i][j],
                                   CMFutureCorrectedMeanRad[i][j], CMFutureCorrectedMeanWindS[i][j],
                                   CMFutureCorrectedStdTdryb[i][j], CMFutureCorrectedStdPressure[i][j],
                                   CMFutureCorrectedStdRad[i][j], CMFutureCorrectedStdWindS[i][j]))

FileMeanStdCMFutureCorrection.close()

print('\nFinished calculation of monthly mean/std of variables for bias corrected future CM data.')
print('Wrote data to file: ', FileNameMeanStdCMFutureCorrection)

# Calculate mean and standard deviation of daily weather variables from EPW for each month of the past period
EPWPastMeanTdryb = numpy.zeros((NYearsPast, NMonths))
EPWPastMeanPressure = numpy.zeros((NYearsPast, NMonths))
EPWPastMeanRad = numpy.zeros((NYearsPast, NMonths))
EPWPastMeanWindS = numpy.zeros((NYearsPast, NMonths))

EPWPastStdTdryb = numpy.zeros((NYearsPast, NMonths))
EPWPastStdPressure = numpy.zeros((NYearsPast, NMonths))
EPWPastStdRad = numpy.zeros((NYearsPast, NMonths))
EPWPastStdWindS = numpy.zeros((NYearsPast, NMonths))

for i in range(0, NYearsPast):
    for j in range(0, NMonths):
        EPWPastMeanTdryb[i][j] = numpy.sum(EPWPastTdrybQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        EPWPastMeanPressure[i][j] = numpy.sum(EPWPastPressureQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        EPWPastMeanRad[i][j] = numpy.sum(EPWPastRadQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]
        EPWPastMeanWindS[i][j] = numpy.sum(EPWPastWindSQuantiles[i][j][0:QuantileSteps[j]])/QuantileSteps[j]

        EPWPastStdTdryb[i][j] = numpy.sqrt(numpy.sum(
            (EPWPastTdrybQuantiles[i][j][0:QuantileSteps[j]] - EPWPastMeanTdryb[i][j]) ** 2) /
                                                    QuantileSteps[j])
        EPWPastStdPressure[i][j] = numpy.sqrt(numpy.sum(
            (EPWPastPressureQuantiles[i][j][0:QuantileSteps[j]] - EPWPastMeanPressure[i][j]) ** 2) /
                                           QuantileSteps[j])
        EPWPastStdRad[i][j] = numpy.sqrt(numpy.sum(
            (EPWPastRadQuantiles[i][j][0:QuantileSteps[j]] - EPWPastMeanRad[i][j]) ** 2) /
                                           QuantileSteps[j])
        EPWPastStdWindS[i][j] = numpy.sqrt(numpy.sum(
            (EPWPastWindSQuantiles[i][j][0:QuantileSteps[j]] - EPWPastMeanWindS[i][j]) ** 2) /
                                           QuantileSteps[j])


# Save data to file
FileMeanStdEPWPast = open(FileNameMeanStdEPWPast, "w")
FileMeanStdEPWPast.write("# 0: Year [-] \t 1: Month [-] \t 2: MeanTdryb [K] \t 3: MeanPressure [Pa] \t 4: MeanRad [W m-2] \t 5: MeanWindS [m s-1]"
                                   "\t 6: StdTdryb [K] \t 7: StdPressure [Pa] \t 8: StdRad [W m-2] \t 9: StdWindS [m s-1] \n")

for i in range(0, NYearsPast):
    for j in range(0, NMonths):
        FileMeanStdEPWPast.write("%i \t %i \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n"
                                % (i, j, EPWPastMeanTdryb[i][j], EPWPastMeanPressure[i][j],
                                   EPWPastMeanRad[i][j], EPWPastMeanWindS[i][j],
                                   EPWPastStdTdryb[i][j], EPWPastStdPressure[i][j],
                                   EPWPastStdRad[i][j], EPWPastStdWindS[i][j]))

FileMeanStdEPWPast.close()

print('\nFinished calculation of monthly mean/std of variables for past EPW data.')
print('Wrote data to file: ', FileNameMeanStdEPWPast)
