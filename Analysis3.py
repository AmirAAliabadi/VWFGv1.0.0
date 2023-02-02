# Quantile-Quantile analysis of CM versus EPW data
# Output matching indices for year and month of past EPW data corresponding to CM for validation period
# Output matching indices for year and month of past EPW data corresponding to CM for future period
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
PathEPWHistorical = ("C:/GoogleDrive/U Guelph/Projects/VWFG/EPW/EPW-Historical")
PathEPWValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/EPW/EPW-Validation")
PathCMValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/CM-Validation-Bias-Corrected")
PathCMFuture = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/CM-Future-Bias-Corrected")
PathResults = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/")

FileNameMatchCMEPWValid = PathResults + 'Match-CM-EPW-Validation.txt'
FileNameMatchCMEPWFuture = PathResults + 'Match-CM-EPW-Future.txt'

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

# Loading and analyzing EPW and CM Data
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

# Read bias corrected climate data for the future period

CMFutureTdrybQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFuturePressureQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFutureRadQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))
CMFutureWindSQuantiles = numpy.zeros((NYearsFuture, NMonths, QuantileStepsMax))

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

    for l in range(0, DailyNPoints):
        CMFutureTdrybQuantiles[i][int(data[l][0])][int(data[l][1])] = data[l][2]
        CMFuturePressureQuantiles[i][int(data[l][0])][int(data[l][1])] = data[l][3]
        CMFutureRadQuantiles[i][int(data[l][0])][int(data[l][1])] = data[l][4]
        CMFutureWindSQuantiles[i][int(data[l][0])][int(data[l][1])] = data[l][5]
    i = i + 1

print('\nFinished reading quantile-quantile future CM data.')
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

# Calculate Finkelstein-Schafer (FS) Method for future period
# Define matrix for FS for each validation year and month
FSFutureTdryb = numpy.zeros((NYearsFuture, NMonths, NYearsPast, NMonths))
FSFuturePressure = numpy.zeros((NYearsFuture, NMonths, NYearsPast, NMonths))
FSFutureRad = numpy.zeros((NYearsFuture, NMonths, NYearsPast, NMonths))
FSFutureWindS = numpy.zeros((NYearsFuture, NMonths, NYearsPast, NMonths))
WSFuture = numpy.zeros((NYearsFuture, NMonths, NYearsPast, NMonths))

for i in range(0, NYearsFuture):
        for j in range(0, NMonths):
            for k in range(0, NYearsPast):
                for l in range(0, NMonths):
                    FSFutureTdryb[i][j][k][l] = numpy.nansum(
                        numpy.abs(CMFutureTdrybQuantiles[i][j][:] - EPWPastTdrybQuantiles[k][l][:])) / QuantileSteps[j]
                    FSFuturePressure[i][j][k][l] = numpy.nansum(
                        numpy.abs(CMFuturePressureQuantiles[i][j][:] - EPWPastPressureQuantiles[k][l][:])) / QuantileSteps[j]
                    FSFutureRad[i][j][k][l] = numpy.nansum(
                        numpy.abs(CMFutureRadQuantiles[i][j][:] - EPWPastRadQuantiles[k][l][:])) / QuantileSteps[j]
                    FSFutureWindS[i][j][k][l] = numpy.nansum(
                        numpy.abs(CMFutureWindSQuantiles[i][j][:] - EPWPastWindSQuantiles[k][l][:])) / QuantileSteps[j]
                    WSFuture[i][j][k][l] = FSWeightsTdryb[j] * FSFutureTdryb[i][j][k][l] + FSWeightsPressure[j] * FSFuturePressure[i][j][k][l] + \
                       FSWeightsRad[j] * FSFutureRad[i][j][k][l] + FSWeightsWindS[j] * FSFutureWindS[i][j][k][l]

print('\nApplied Finkelstein-Schafer (FS) method for validation and future CM data.')

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

# Look for minimum weighted sum for every month and year of the future period
# Check within all the past months to find the lowest weighted sum
FutureIndex = numpy.zeros((NYearsFuture, NMonths, 2))

for i in range(0, NYearsFuture):
        for j in range(0, NMonths):
            MinimumWS = 10e6
            for k in range(0, NYearsPast):
                for l in range(0, NMonths):
                    if WSFuture[i][j][k][l] < MinimumWS:
                        MinimumWS = WSFuture[i][j][k][l]
                        FutureIndex[i][j][0] = k
                        FutureIndex[i][j][1] = l

print('\nCalculated year/month indices for matching EPW records for validation and future CM data.')

# Write matching years and months to file

# Change the directory
os.chdir(PathResults)

FileMatchCMEPWValid = open(FileNameMatchCMEPWValid, "w")
FileMatchCMEPWValid.write("# 0: Valid Year [-] \t 1: Valid Month [-] \t 2: Match Year [-] \t 3: Match Month [-] \n")

for i in range(0, NYearsValid):
        for j in range(0, NMonths):
            FileMatchCMEPWValid.write("%i \t %i \t %i \t %i \n" % (i, j, ValidIndex[i][j][0], ValidIndex[i][j][1]))

FileMatchCMEPWValid.close()

print('\nWrote year/month indices for matching EPW records for validation CM data.')
print('Wrote data to file: ', FileNameMatchCMEPWValid)

FileMatchCMEPWFuture = open(FileNameMatchCMEPWFuture, "w")
FileMatchCMEPWFuture.write("# 0: Future Year [-] \t 1: Future Month [-] \t 2: Match Year [-] \t 3: Match Month [-] \n")

for i in range(0, NYearsFuture):
        for j in range(0, NMonths):
            FileMatchCMEPWFuture.write("%i \t %i \t %i \t %i \n" % (i, j, FutureIndex[i][j][0], FutureIndex[i][j][1]))

FileMatchCMEPWFuture.close()

print('\nWrote year/month indices for matching EPW records for future CM data.')
print('Wrote data to file: ', FileNameMatchCMEPWFuture)