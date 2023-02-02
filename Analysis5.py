# Shift/stretch weather variables for the morphed EPW files for validation and future periods
# Computing errors for the method before and after shifting/stretching for the validation period
# Updating the morphed EPW files by shifting and stretching
# Atmospheric Innovations Research (AIR) Laboratory, University of Guelph
# Developed by Amir A. Aliabadi and Rachel M. McLeod
# Last updated: 2022-04-05

import numpy
from array import array
import math
import pandas
import os
import urllib.request
import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import shutil

# Constants of simulation
########################################################################################################################

# Historical period: 1980-1999
# Validation period: 2007-2020
# Future period: 2021-2100

# Define directories and file names
PathMain = ("C:/GoogleDrive/U Guelph/Projects/VWFG")
PathEPWHistorical = ("C:/GoogleDrive/U Guelph/Projects/VWFG/EPW/EPW-Historical")
PathEPWValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/EPW/EPW-Validation")
PathCMValidationMatch = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results")
PathCMValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/CM-Validation-Bias-Corrected")
PathCMFuture = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/CM-Future-Bias-Corrected")
PathMatchValidation = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/Match-CM-EPW-Validation")
PathEPWFutureNotCorrected = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/EPW-Future-Not-Corrected")
PathEPWFutureCorrected = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/EPW-Future-Corrected")
PathResults = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/")

FileNameErrorsHourlyMultiYearValidNoShiftNoStretch = PathResults + 'Errors-Morphed-Valid-No-Shift-No-Stretch.txt'
FileNameErrorsHourlyMultiYearValidShiftedStretched = PathResults + 'Errors-Morphed-Valid-Shifted-Stretched.txt'

MatchValidFile = "Match-CM-EPW-Validation.txt"
MatchFutureFile = "Match-CM-EPW-Future.txt"
MeanStdEPWPastFile = "Mean-Std-EPW-Past.txt"
MeanStdCMValidFile = "Mean-Std-CM-Validation-Correction.txt"
MeanStdCMFutureFile = "Mean-Std-CM-Future-Correction.txt"
EPWMorphedShiftedStretchedFile = "EPW-Morphed-Shifted-Stretched-Year-"

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
# Soil layers information in EPW file
SoilLayersN = 3
SoilTemperatureLineN = 4
SoilTemperatureStart = 6
SoilTemperatureJump = 16

# Define minimum wind speed [m s-1]
MinimumWindS = 0.1

# Calculations
########################################################################################################################
# Loading and analyzing EPW data for the validation period
# Go to results directory
os.chdir(PathResults)

# Load matching indices for the validation period
MatchingIndicesValid = numpy.loadtxt(MatchValidFile)

# Go to directory that contains the historic EPW file
os.chdir(PathEPWHistorical)
# Remember the file names in order in this directory
EPWHistoricalFiles = os.listdir(PathEPWHistorical)

# Read matching EPW hourly data and store values to a matrix
MatchingEPWValidTdryb = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
MatchingEPWValidPressure = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
MatchingEPWValidRad = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
MatchingEPWValidWindS = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))

for i in MatchingIndicesValid[:, :]:
    ValidationYear = int(i[0])                   # Current validation year
    ValidationMonth = int(i[1])                  # Current validation month
    HistoricalYear = int(i[2])                   # Matching year in historical record
    StartRow = MonthIndexEPWPandas[int(i[3])]          # Defines what rows we should start reading from the EPW file
    EndRow = MonthIndexEPWPandas[int(i[3] + 1)]        # Define what row we should be ending at reading from the EPW file
    YearFile = EPWHistoricalFiles[HistoricalYear]

    data = pandas.read_csv(YearFile, header = NSkipHeader)
    data = data.values

    # Note pandas.read_csv() cannot read first hour of data in EPW file
    # Consider the second hour data as first hour data
    if StartRow == 0:
        MatchingEPWValidTdryb[ValidationYear][ValidationMonth][0] = CelciusToKelvin + data[0, 6]
        MatchingEPWValidPressure[ValidationYear][ValidationMonth][0] = data[0, 9]
        # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
        MatchingEPWValidRad[ValidationYear][ValidationMonth][0] = data[0, 12] + data[0, 14] + data[0, 15]
        MatchingEPWValidWindS[ValidationYear][ValidationMonth][0] = data[0, 21]
        Index = 1
    # Otherwise fine
    else:
        Index = 0

    for l in range(StartRow, EndRow):
        MatchingEPWValidTdryb[ValidationYear][ValidationMonth][Index] = CelciusToKelvin + data[l, 6]
        MatchingEPWValidPressure[ValidationYear][ValidationMonth][Index] = data[l, 9]
        # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
        MatchingEPWValidRad[ValidationYear][ValidationMonth][Index] = data[l, 12] + data[l, 14] + data[l, 15]
        MatchingEPWValidWindS[ValidationYear][ValidationMonth][Index] = data[l, 21]
        Index = Index + 1

print('\nFinished reading matching historical EPW data for validation period.')

# Go to directory that contains the validation EPW file
os.chdir(PathEPWValidation)
# Remember the file names in order in this directory
EPWValidationFiles = os.listdir(PathEPWValidation)

# Read matching EPW hourly data and store values to a matrix
EPWValidTdryb = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
EPWValidPressure = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
EPWValidRad = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
EPWValidWindS = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))

for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        YearFile = EPWValidationFiles[i]
        data = pandas.read_csv(YearFile, header = NSkipHeader)
        data = data.values

        # Note pandas.read_csv() cannot read first hour of data in EPW file
        # Consider the second hour data as first hour data
        if MonthIndexEPWPandas[j] == 0:
            EPWValidTdryb[i][j][0] = CelciusToKelvin + data[0, 6]
            EPWValidPressure[i][j][0] = data[0, 9]
            # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
            EPWValidRad[i][j][0] = data[0, 12] + data[0, 14] + data[0, 15]
            EPWValidWindS[i][j][0] = data[0, 21]
            Index = 1
        # Otherwise fine
        else:
            Index = 0

        for l in range(MonthIndexEPWPandas[j], MonthIndexEPWPandas[j+1]):
            EPWValidTdryb[i][j][Index] = CelciusToKelvin + data[l, 6]
            EPWValidPressure[i][j][Index] = data[l, 9]
            # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
            EPWValidRad[i][j][Index] = data[l, 12] + data[l, 14] + data[l, 15]
            EPWValidWindS[i][j][Index] = data[l, 21]
            Index = Index + 1

print('\nFinished reading validation EPW data.')

# Calculate month by month bias and RMSE between matched EPW dataset and the actual EPW dataset for the validation period

BiasHourlyValidTdryb = numpy.zeros((NYearsValid, NMonths))
BiasHourlyValidPressure = numpy.zeros((NYearsValid, NMonths))
BiasHourlyValidRad = numpy.zeros((NYearsValid, NMonths))
BiasHourlyValidWindS = numpy.zeros((NYearsValid, NMonths))

RMSEHourlyValidTdryb = numpy.zeros((NYearsValid, NMonths))
RMSEHourlyValidPressure = numpy.zeros((NYearsValid, NMonths))
RMSEHourlyValidRad = numpy.zeros((NYearsValid, NMonths))
RMSEHourlyValidWindS = numpy.zeros((NYearsValid, NMonths))

BiasHourlyMultiYearValidTdryb = numpy.zeros((NMonths,1))
BiasHourlyMultiYearValidPressure = numpy.zeros((NMonths, 1))
BiasHourlyMultiYearValidRad = numpy.zeros((NMonths, 1))
BiasHourlyMultiYearValidWindS = numpy.zeros((NMonths, 1))

RMSEHourlyMultiYearValidTdryb = numpy.zeros((NMonths, 1))
RMSEHourlyMultiYearValidPressure = numpy.zeros((NMonths, 1))
RMSEHourlyMultiYearValidRad = numpy.zeros((NMonths, 1))
RMSEHourlyMultiYearValidWindS = numpy.zeros((NMonths, 1))

for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        BiasHourlyValidTdryb[i][j] = numpy.sum(MatchingEPWValidTdryb[i][j][0:MonthIndex[j]] - EPWValidTdryb[i][j][0:MonthIndex[j]])/(MonthIndex[j])
        RMSEHourlyValidTdryb[i][j] = numpy.sqrt(numpy.sum((MatchingEPWValidTdryb[i][j][0:MonthIndex[j]] - EPWValidTdryb[i][j][0:MonthIndex[j]])**2)/(MonthIndex[j]))
        BiasHourlyValidPressure[i][j] = numpy.sum(MatchingEPWValidPressure[i][j][0:MonthIndex[j]] - EPWValidPressure[i][j][0:MonthIndex[j]]) /(MonthIndex[j])
        RMSEHourlyValidPressure[i][j] = numpy.sqrt(numpy.sum((MatchingEPWValidPressure[i][j][0:MonthIndex[j]] - EPWValidPressure[i][j][0:MonthIndex[j]])**2)/(MonthIndex[j]))
        BiasHourlyValidRad[i][j] = numpy.sum(MatchingEPWValidRad[i][j][0:MonthIndex[j]] - EPWValidRad[i][j][0:MonthIndex[j]]) / (MonthIndex[j])
        RMSEHourlyValidRad[i][j] = numpy.sqrt(numpy.sum((MatchingEPWValidRad[i][j][0:MonthIndex[j]] - EPWValidRad[i][j][0:MonthIndex[j]])**2) / (MonthIndex[j]))
        BiasHourlyValidWindS[i][j] = numpy.sum(MatchingEPWValidWindS[i][j][0:MonthIndex[j]] - EPWValidWindS[i][j][0:MonthIndex[j]]) / (MonthIndex[j])
        RMSEHourlyValidWindS[i][j] = numpy.sqrt(numpy.sum((MatchingEPWValidWindS[i][j][0:MonthIndex[j]] - EPWValidWindS[i][j][0:MonthIndex[j]])**2) / (MonthIndex[j]))

print('\nFinished calculating Bias/RMSE for matching historical vs. validation EPW data.')

# Calculate average Bias and RMSE for each month during the validation period
for j in range(0, NMonths):
    SumBiasTdryb = 0
    SumBiasPressure = 0
    SumBiasRad = 0
    SumBiasWindS = 0
    SumRMSETdryb = 0
    SumRMSEPressure = 0
    SumRMSERad = 0
    SumRMSEWindS = 0
    for i in range(0, NYearsValid):
        SumBiasTdryb = SumBiasTdryb + BiasHourlyValidTdryb[i][j]
        SumBiasPressure = SumBiasPressure + BiasHourlyValidPressure[i][j]
        SumBiasRad = SumBiasRad + BiasHourlyValidRad[i][j]
        SumBiasWindS = SumBiasWindS + BiasHourlyValidWindS[i][j]

        SumRMSETdryb = SumRMSETdryb + RMSEHourlyValidTdryb[i][j]
        SumRMSEPressure = SumRMSEPressure + RMSEHourlyValidPressure[i][j]
        SumRMSERad = SumRMSERad + RMSEHourlyValidRad[i][j]
        SumRMSEWindS = SumRMSEWindS + RMSEHourlyValidWindS[i][j]

    BiasHourlyMultiYearValidTdryb[j] = SumBiasTdryb / NYearsValid
    BiasHourlyMultiYearValidPressure[j] = SumBiasPressure / NYearsValid
    BiasHourlyMultiYearValidRad[j] = SumBiasRad / NYearsValid
    BiasHourlyMultiYearValidWindS[j] = SumBiasWindS / NYearsValid

    RMSEHourlyMultiYearValidTdryb[j] = SumRMSETdryb / NYearsValid
    RMSEHourlyMultiYearValidPressure[j] = SumRMSEPressure / NYearsValid
    RMSEHourlyMultiYearValidRad[j] = SumRMSERad / NYearsValid
    RMSEHourlyMultiYearValidWindS[j] = SumRMSEWindS / NYearsValid

print('Finished calculating multi-year Bias/RMSE for matching historical vs. validation EPW data.')

# Save to file
FileErrorsHourlyMultiYearValidNoShiftNoStretch = open(FileNameErrorsHourlyMultiYearValidNoShiftNoStretch, "w")
FileErrorsHourlyMultiYearValidNoShiftNoStretch.write("# 0: Month [-] \t 1: BiasTdryb [K] \t 2: BiasPressure [Pa] \t 3: BiasRad [W m-2] \t 4: BiasWindS [m s-1] \t"
                                   "5: RMSETdryb [K] \t 6: RMSEPressure [Pa] \t 7: RMSERad [W m-2] \t 8: RMSEWindS [m s-1] \n")

for j in range(0, NMonths):
    FileErrorsHourlyMultiYearValidNoShiftNoStretch.write("%i \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n"
                             % (j, BiasHourlyMultiYearValidTdryb[j], BiasHourlyMultiYearValidPressure[j],
                                BiasHourlyMultiYearValidRad[j], BiasHourlyMultiYearValidWindS[j],
                                RMSEHourlyMultiYearValidTdryb[j], RMSEHourlyMultiYearValidPressure[j],
                                RMSEHourlyMultiYearValidRad[j], RMSEHourlyMultiYearValidWindS[j]))

FileErrorsHourlyMultiYearValidNoShiftNoStretch.close()

print('Wrote data to file: ', FileNameErrorsHourlyMultiYearValidNoShiftNoStretch)

# Load mean and standard deviation of monthly weather variables from CM dataset for the validation period
# Go to results directory
os.chdir(PathResults)

MeanStdCMValid = numpy.loadtxt(MeanStdCMValidFile)

MeanTdrybCMValid = numpy.zeros((NYearsValid, NMonths))
MeanPressureCMValid = numpy.zeros((NYearsValid, NMonths))
MeanRadCMValid = numpy.zeros((NYearsValid, NMonths))
MeanWindSCMValid = numpy.zeros((NYearsValid, NMonths))
StdTdrybCMValid = numpy.zeros((NYearsValid, NMonths))
StdPressureCMValid = numpy.zeros((NYearsValid, NMonths))
StdRadCMValid = numpy.zeros((NYearsValid, NMonths))
StdWindSCMValid = numpy.zeros((NYearsValid, NMonths))

for i in range(0, int(NYearsValid * NMonths)):
        MeanTdrybCMValid[int(MeanStdCMValid[i][0])][int(MeanStdCMValid[i][1])] = MeanStdCMValid[i][2]
        MeanPressureCMValid[int(MeanStdCMValid[i][0])][int(MeanStdCMValid[i][1])] = MeanStdCMValid[i][3]
        MeanRadCMValid[int(MeanStdCMValid[i][0])][int(MeanStdCMValid[i][1])] = MeanStdCMValid[i][4]
        MeanWindSCMValid[int(MeanStdCMValid[i][0])][int(MeanStdCMValid[i][1])] = MeanStdCMValid[i][5]
        StdTdrybCMValid[int(MeanStdCMValid[i][0])][int(MeanStdCMValid[i][1])] = MeanStdCMValid[i][6]
        StdPressureCMValid[int(MeanStdCMValid[i][0])][int(MeanStdCMValid[i][1])] = MeanStdCMValid[i][7]
        StdRadCMValid[int(MeanStdCMValid[i][0])][int(MeanStdCMValid[i][1])] = MeanStdCMValid[i][8]
        StdWindSCMValid[int(MeanStdCMValid[i][0])][int(MeanStdCMValid[i][1])] = MeanStdCMValid[i][9]

# Load mean and standard deviation of monthly weather variables for EPW dataset during the past period
MeanStdEPWPast = numpy.loadtxt(MeanStdEPWPastFile)

MeanTdrybEPWPast = numpy.zeros((NYearsPast, NMonths))
MeanPressureEPWPast = numpy.zeros((NYearsPast, NMonths))
MeanRadEPWPast = numpy.zeros((NYearsPast, NMonths))
MeanWindSEPWPast = numpy.zeros((NYearsPast, NMonths))
StdTdrybEPWPast = numpy.zeros((NYearsPast, NMonths))
StdPressureEPWPast = numpy.zeros((NYearsPast, NMonths))
StdRadEPWPast = numpy.zeros((NYearsPast, NMonths))
StdWindSEPWPast = numpy.zeros((NYearsPast, NMonths))

for i in range(0, int(NYearsPast * NMonths)):
        MeanTdrybEPWPast[int(MeanStdEPWPast[i][0])][int(MeanStdEPWPast[i][1])] = MeanStdEPWPast[i][2]
        MeanPressureEPWPast[int(MeanStdEPWPast[i][0])][int(MeanStdEPWPast[i][1])] = MeanStdEPWPast[i][3]
        MeanRadEPWPast[int(MeanStdEPWPast[i][0])][int(MeanStdEPWPast[i][1])] = MeanStdEPWPast[i][4]
        MeanWindSEPWPast[int(MeanStdEPWPast[i][0])][int(MeanStdEPWPast[i][1])] = MeanStdEPWPast[i][5]
        StdTdrybEPWPast[int(MeanStdEPWPast[i][0])][int(MeanStdEPWPast[i][1])] = MeanStdEPWPast[i][6]
        StdPressureEPWPast[int(MeanStdEPWPast[i][0])][int(MeanStdEPWPast[i][1])] = MeanStdEPWPast[i][7]
        StdRadEPWPast[int(MeanStdEPWPast[i][0])][int(MeanStdEPWPast[i][1])] = MeanStdEPWPast[i][8]
        StdWindSEPWPast[int(MeanStdEPWPast[i][0])][int(MeanStdEPWPast[i][1])] = MeanStdEPWPast[i][9]

# Calculate monthly mean and standard deviation of weather variables for each of the matched EPW files from the historical period for the validation period
MatchingEPWValidMeanTdryb = numpy.zeros((NYearsValid, NMonths))
MatchingEPWValidMeanPressure = numpy.zeros((NYearsValid, NMonths))
MatchingEPWValidMeanRad = numpy.zeros((NYearsValid, NMonths))
MatchingEPWValidMeanWindS = numpy.zeros((NYearsValid, NMonths))

MatchingEPWValidStdTdryb = numpy.zeros((NYearsValid, NMonths))
MatchingEPWValidStdPressure = numpy.zeros((NYearsValid, NMonths))
MatchingEPWValidStdRad = numpy.zeros((NYearsValid, NMonths))
MatchingEPWValidStdWindS = numpy.zeros((NYearsValid, NMonths))

for i in range(0, int(NYearsValid * NMonths)):
    MatchingEPWValidMeanTdryb[int(MatchingIndicesValid[i][0])][int(MatchingIndicesValid[i][1])] = MeanTdrybEPWPast[int(MatchingIndicesValid[i][2])][int(MatchingIndicesValid[i][3])]
    MatchingEPWValidMeanPressure[int(MatchingIndicesValid[i][0])][int(MatchingIndicesValid[i][1])] = MeanPressureEPWPast[int(MatchingIndicesValid[i][2])][int(MatchingIndicesValid[i][3])]
    MatchingEPWValidMeanRad[int(MatchingIndicesValid[i][0])][int(MatchingIndicesValid[i][1])] = MeanRadEPWPast[int(MatchingIndicesValid[i][2])][int(MatchingIndicesValid[i][3])]
    MatchingEPWValidMeanWindS[int(MatchingIndicesValid[i][0])][int(MatchingIndicesValid[i][1])] = MeanWindSEPWPast[int(MatchingIndicesValid[i][2])][int(MatchingIndicesValid[i][3])]

    MatchingEPWValidStdTdryb[int(MatchingIndicesValid[i][0])][int(MatchingIndicesValid[i][1])] = StdTdrybEPWPast[int(MatchingIndicesValid[i][2])][int(MatchingIndicesValid[i][3])]
    MatchingEPWValidStdPressure[int(MatchingIndicesValid[i][0])][int(MatchingIndicesValid[i][1])] = StdPressureEPWPast[int(MatchingIndicesValid[i][2])][int(MatchingIndicesValid[i][3])]
    MatchingEPWValidStdRad[int(MatchingIndicesValid[i][0])][int(MatchingIndicesValid[i][1])] = StdRadEPWPast[int(MatchingIndicesValid[i][2])][int(MatchingIndicesValid[i][3])]
    MatchingEPWValidStdWindS[int(MatchingIndicesValid[i][0])][int(MatchingIndicesValid[i][1])] = StdWindSEPWPast[int(MatchingIndicesValid[i][2])][int(MatchingIndicesValid[i][3])]

# Use the matching EPW and CM model means and standard deviations to stretch and shift the data for the validation period
MatchingEPWValidTdrybShiftedStretched = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
MatchingEPWValidPressureShiftedStretched = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
MatchingEPWValidRadShiftedStretched = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))
MatchingEPWValidWindSShiftedStretched = numpy.zeros((NYearsValid, NMonths, NHoursInMonth))

for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        MatchingEPWValidTdrybShiftedStretched[i][j][:] = MatchingEPWValidMeanTdryb[i][j] + \
            MeanTdrybCMValid[i][j] - MatchingEPWValidMeanTdryb[i][j] + \
            (MatchingEPWValidTdryb[i][j][:] - MatchingEPWValidMeanTdryb[i][j]) * \
            StdTdrybCMValid[i][j] / MatchingEPWValidStdTdryb[i][j]
        MatchingEPWValidPressureShiftedStretched[i][j][:] = MatchingEPWValidMeanPressure[i][j] + \
            MeanPressureCMValid[i][j] - MatchingEPWValidMeanPressure[i][j] + \
            (MatchingEPWValidPressure[i][j][:] - MatchingEPWValidMeanPressure[i][j]) * \
            StdPressureCMValid[i][j] / MatchingEPWValidStdPressure[i][j]
        MatchingEPWValidRadShiftedStretched[i][j][:] = MatchingEPWValidMeanRad[i][j] + \
            MeanRadCMValid[i][j] - MatchingEPWValidMeanRad[i][j] + \
            (MatchingEPWValidRad[i][j][:] - MatchingEPWValidMeanRad[i][j]) * \
            StdRadCMValid[i][j] / MatchingEPWValidStdRad[i][j]
        MatchingEPWValidWindSShiftedStretched[i][j][:] = MatchingEPWValidMeanWindS[i][j] + \
            MeanWindSCMValid[i][j] - MatchingEPWValidMeanWindS[i][j] + \
            (MatchingEPWValidWindS[i][j][:] - MatchingEPWValidMeanWindS[i][j]) * \
            StdWindSCMValid[i][j] / MatchingEPWValidStdWindS[i][j]

# Check for wind speed not to be negative, otherwise apply a minimum wind speed
for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        for k in range(0, NHoursInMonth):
            MatchingEPWValidWindSShiftedStretched[i][j][k] = max(MinimumWindS, MatchingEPWValidWindSShiftedStretched[i][j][k])

# Calculate month by month bias and RMSE between matched, shifted, and stretched EPW dataset and the actual EPW dataset for the validation period

BiasHourlyValidTdrybShiftedStretched = numpy.zeros((NYearsValid, NMonths))
BiasHourlyValidPressureShiftedStretched = numpy.zeros((NYearsValid, NMonths))
BiasHourlyValidRadShiftedStretched = numpy.zeros((NYearsValid, NMonths))
BiasHourlyValidWindSShiftedStretched = numpy.zeros((NYearsValid, NMonths))

RMSEHourlyValidTdrybShiftedStretched = numpy.zeros((NYearsValid, NMonths))
RMSEHourlyValidPressureShiftedStretched = numpy.zeros((NYearsValid, NMonths))
RMSEHourlyValidRadShiftedStretched = numpy.zeros((NYearsValid, NMonths))
RMSEHourlyValidWindSShiftedStretched = numpy.zeros((NYearsValid, NMonths))

BiasHourlyMultiYearValidTdrybShiftedStretched = numpy.zeros((NMonths,1))
BiasHourlyMultiYearValidPressureShiftedStretched = numpy.zeros((NMonths, 1))
BiasHourlyMultiYearValidRadShiftedStretched = numpy.zeros((NMonths, 1))
BiasHourlyMultiYearValidWindSShiftedStretched = numpy.zeros((NMonths, 1))

RMSEHourlyMultiYearValidTdrybShiftedStretched = numpy.zeros((NMonths, 1))
RMSEHourlyMultiYearValidPressureShiftedStretched = numpy.zeros((NMonths, 1))
RMSEHourlyMultiYearValidRadShiftedStretched = numpy.zeros((NMonths, 1))
RMSEHourlyMultiYearValidWindSShiftedStretched = numpy.zeros((NMonths, 1))

for i in range(0, NYearsValid):
    for j in range(0, NMonths):
        BiasHourlyValidTdrybShiftedStretched[i][j] = numpy.sum(MatchingEPWValidTdrybShiftedStretched[i][j][0:MonthIndex[j]] - EPWValidTdryb[i][j][0:MonthIndex[j]])/(MonthIndex[j])
        RMSEHourlyValidTdrybShiftedStretched[i][j] = numpy.sqrt(numpy.sum((MatchingEPWValidTdrybShiftedStretched[i][j][0:MonthIndex[j]] - EPWValidTdryb[i][j][0:MonthIndex[j]])**2)/(MonthIndex[j]))
        BiasHourlyValidPressureShiftedStretched[i][j] = numpy.sum(MatchingEPWValidPressureShiftedStretched[i][j][0:MonthIndex[j]] - EPWValidPressure[i][j][0:MonthIndex[j]]) /(MonthIndex[j])
        RMSEHourlyValidPressureShiftedStretched[i][j] = numpy.sqrt(numpy.sum((MatchingEPWValidPressureShiftedStretched[i][j][0:MonthIndex[j]] - EPWValidPressure[i][j][0:MonthIndex[j]])**2)/(MonthIndex[j]))
        BiasHourlyValidRadShiftedStretched[i][j] = numpy.sum(MatchingEPWValidRadShiftedStretched[i][j][0:MonthIndex[j]] - EPWValidRad[i][j][0:MonthIndex[j]]) / (MonthIndex[j])
        RMSEHourlyValidRadShiftedStretched[i][j] = numpy.sqrt(numpy.sum((MatchingEPWValidRadShiftedStretched[i][j][0:MonthIndex[j]] - EPWValidRad[i][j][0:MonthIndex[j]])**2) / (MonthIndex[j]))
        BiasHourlyValidWindSShiftedStretched[i][j] = numpy.sum(MatchingEPWValidWindSShiftedStretched[i][j][0:MonthIndex[j]] - EPWValidWindS[i][j][0:MonthIndex[j]]) / (MonthIndex[j])
        RMSEHourlyValidWindSShiftedStretched[i][j] = numpy.sqrt(numpy.sum((MatchingEPWValidWindSShiftedStretched[i][j][0:MonthIndex[j]] - EPWValidWindS[i][j][0:MonthIndex[j]])**2) / (MonthIndex[j]))

print('\nFinished calculating Bias/RMSE for corrected matching historical vs. validation EPW data.')

# Calculate average Bias and RMSE for each month during the validation period
for j in range(0, NMonths):
    SumBiasTdryb = 0
    SumBiasPressure = 0
    SumBiasRad = 0
    SumBiasWindS = 0
    SumRMSETdryb = 0
    SumRMSEPressure = 0
    SumRMSERad = 0
    SumRMSEWindS = 0
    for i in range(0, NYearsValid):
        SumBiasTdryb = SumBiasTdryb + BiasHourlyValidTdrybShiftedStretched[i][j]
        SumBiasPressure = SumBiasPressure + BiasHourlyValidPressureShiftedStretched[i][j]
        SumBiasRad = SumBiasRad + BiasHourlyValidRadShiftedStretched[i][j]
        SumBiasWindS = SumBiasWindS + BiasHourlyValidWindSShiftedStretched[i][j]

        SumRMSETdryb = SumRMSETdryb + RMSEHourlyValidTdrybShiftedStretched[i][j]
        SumRMSEPressure = SumRMSEPressure + RMSEHourlyValidPressureShiftedStretched[i][j]
        SumRMSERad = SumRMSERad + RMSEHourlyValidRadShiftedStretched[i][j]
        SumRMSEWindS = SumRMSEWindS + RMSEHourlyValidWindSShiftedStretched[i][j]

    BiasHourlyMultiYearValidTdrybShiftedStretched[j] = SumBiasTdryb / NYearsValid
    BiasHourlyMultiYearValidPressureShiftedStretched[j] = SumBiasPressure / NYearsValid
    BiasHourlyMultiYearValidRadShiftedStretched[j] = SumBiasRad / NYearsValid
    BiasHourlyMultiYearValidWindSShiftedStretched[j] = SumBiasWindS / NYearsValid

    RMSEHourlyMultiYearValidTdrybShiftedStretched[j] = SumRMSETdryb / NYearsValid
    RMSEHourlyMultiYearValidPressureShiftedStretched[j] = SumRMSEPressure / NYearsValid
    RMSEHourlyMultiYearValidRadShiftedStretched[j] = SumRMSERad / NYearsValid
    RMSEHourlyMultiYearValidWindSShiftedStretched[j] = SumRMSEWindS / NYearsValid

print('Finished calculating multi-year Bias/RMSE for corrected matching historical vs. validation EPW data.')

# Save to file
FileErrorsHourlyMultiYearValidShiftedStretched = open(FileNameErrorsHourlyMultiYearValidShiftedStretched, "w")
FileErrorsHourlyMultiYearValidShiftedStretched.write("# 0: Month [-] \t 1: BiasTdryb [K] \t 2: BiasPressure [Pa] \t 3: BiasRad [W m-2] \t 4: BiasWindS [m s-1] \t"
                                   "5: RMSETdryb [K] \t 6: RMSEPressure [Pa] \t 7: RMSERad [W m-2] \t 8: RMSEWindS [m s-1] \n")

for j in range(0, NMonths):
    FileErrorsHourlyMultiYearValidShiftedStretched.write("%i \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n"
                             % (j, BiasHourlyMultiYearValidTdrybShiftedStretched[j], BiasHourlyMultiYearValidPressureShiftedStretched[j],
                                BiasHourlyMultiYearValidRadShiftedStretched[j], BiasHourlyMultiYearValidWindSShiftedStretched[j],
                                RMSEHourlyMultiYearValidTdrybShiftedStretched[j], RMSEHourlyMultiYearValidPressureShiftedStretched[j],
                                RMSEHourlyMultiYearValidRadShiftedStretched[j], RMSEHourlyMultiYearValidWindSShiftedStretched[j]))

FileErrorsHourlyMultiYearValidShiftedStretched.close()

print('Wrote data to file: ', FileNameErrorsHourlyMultiYearValidShiftedStretched)

# Loading and analyzing EPW data for the future period
# Go to results directory
os.chdir(PathResults)

# Load matching indices for the future period
MatchingIndicesFuture = numpy.loadtxt(MatchFutureFile)

# Go to directory that contains the historic EPW file
os.chdir(PathEPWHistorical)
# Remember the file names in order in this directory
EPWHistoricalFiles = os.listdir(PathEPWHistorical)

# Read matching EPW hourly data and store values to a matrix
MatchingEPWFutureTdryb = numpy.zeros((NYearsFuture, NMonths, NHoursInMonth))
MatchingEPWFuturePressure = numpy.zeros((NYearsFuture, NMonths, NHoursInMonth))
MatchingEPWFutureRad = numpy.zeros((NYearsFuture, NMonths, NHoursInMonth))
MatchingEPWFutureWindS = numpy.zeros((NYearsFuture, NMonths, NHoursInMonth))

for i in MatchingIndicesFuture[:, :]:
    FutureYear = int(i[0])                       # Current future year
    FutureMonth = int(i[1])                      # Current future month
    HistoricalYear = int(i[2])                   # Matching year in historical record
    StartRow = MonthIndexEPWPandas[int(i[3])]          # Defines what rows we should start reading from the EPW file
    EndRow = MonthIndexEPWPandas[int(i[3] + 1)]        # Define what row we should be ending at reading from the EPW file
    YearFile = EPWHistoricalFiles[HistoricalYear]

    data = pandas.read_csv(YearFile, header = NSkipHeader)
    data = data.values

    # Note pandas.read_csv() cannot read first hour of data in EPW file
    # Consider the second hour data as first hour data
    if StartRow == 0:
        MatchingEPWFutureTdryb[FutureYear][FutureMonth][0] = CelciusToKelvin + data[0, 6]
        MatchingEPWFuturePressure[FutureYear][FutureMonth][0] = data[0, 9]
        # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
        MatchingEPWFutureRad[FutureYear][FutureMonth][0] = data[0, 12] + data[0, 14] + data[0, 15]
        MatchingEPWFutureWindS[FutureYear][FutureMonth][0] = data[0, 21]
        Index = 1
    # Otherwise fine
    else:
        Index = 0

    for l in range(StartRow, EndRow):
        MatchingEPWFutureTdryb[FutureYear][FutureMonth][Index] = CelciusToKelvin + data[l, 6]
        MatchingEPWFuturePressure[FutureYear][FutureMonth][Index] = data[l, 9]
        # For total global horizontal radiation flux add longwave (12), direct shortwave (14), and diffuse shortwave (15) fluxes
        MatchingEPWFutureRad[FutureYear][FutureMonth][Index] = data[l, 12] + data[l, 14] + data[l, 15]
        MatchingEPWFutureWindS[FutureYear][FutureMonth][Index] = data[l, 21]
        Index = Index + 1

print('\nFinished reading matching historical EPW data for future period.')

# Load mean and standard deviation of monthly weather variables from CM dataset for the future period
# Go to results directory
os.chdir(PathResults)

MeanStdCMFuture = numpy.loadtxt(MeanStdCMFutureFile)

MeanTdrybCMFuture = numpy.zeros((NYearsFuture, NMonths))
MeanPressureCMFuture = numpy.zeros((NYearsFuture, NMonths))
MeanRadCMFuture = numpy.zeros((NYearsFuture, NMonths))
MeanWindSCMFuture = numpy.zeros((NYearsFuture, NMonths))
StdTdrybCMFuture = numpy.zeros((NYearsFuture, NMonths))
StdPressureCMFuture = numpy.zeros((NYearsFuture, NMonths))
StdRadCMFuture = numpy.zeros((NYearsFuture, NMonths))
StdWindSCMFuture = numpy.zeros((NYearsFuture, NMonths))

for i in range(0, int(NYearsFuture * NMonths)):
        MeanTdrybCMFuture[int(MeanStdCMFuture[i][0])][int(MeanStdCMFuture[i][1])] = MeanStdCMFuture[i][2]
        MeanPressureCMFuture[int(MeanStdCMFuture[i][0])][int(MeanStdCMFuture[i][1])] = MeanStdCMFuture[i][3]
        MeanRadCMFuture[int(MeanStdCMFuture[i][0])][int(MeanStdCMFuture[i][1])] = MeanStdCMFuture[i][4]
        MeanWindSCMFuture[int(MeanStdCMFuture[i][0])][int(MeanStdCMFuture[i][1])] = MeanStdCMFuture[i][5]
        StdTdrybCMFuture[int(MeanStdCMFuture[i][0])][int(MeanStdCMFuture[i][1])] = MeanStdCMFuture[i][6]
        StdPressureCMFuture[int(MeanStdCMFuture[i][0])][int(MeanStdCMFuture[i][1])] = MeanStdCMFuture[i][7]
        StdRadCMFuture[int(MeanStdCMFuture[i][0])][int(MeanStdCMFuture[i][1])] = MeanStdCMFuture[i][8]
        StdWindSCMFuture[int(MeanStdCMFuture[i][0])][int(MeanStdCMFuture[i][1])] = MeanStdCMFuture[i][9]

# Calculate monthly mean and standard deviation of weather variables for each of the matched EPW files from the historical period for the future period
MatchingEPWFutureMeanTdryb = numpy.zeros((NYearsFuture, NMonths))
MatchingEPWFutureMeanPressure = numpy.zeros((NYearsFuture, NMonths))
MatchingEPWFutureMeanRad = numpy.zeros((NYearsFuture, NMonths))
MatchingEPWFutureMeanWindS = numpy.zeros((NYearsFuture, NMonths))

MatchingEPWFutureStdTdryb = numpy.zeros((NYearsFuture, NMonths))
MatchingEPWFutureStdPressure = numpy.zeros((NYearsFuture, NMonths))
MatchingEPWFutureStdRad = numpy.zeros((NYearsFuture, NMonths))
MatchingEPWFutureStdWindS = numpy.zeros((NYearsFuture, NMonths))

for i in range(0, int(NYearsFuture * NMonths)):
    MatchingEPWFutureMeanTdryb[int(MatchingIndicesFuture[i][0])][int(MatchingIndicesFuture[i][1])] = MeanTdrybEPWPast[int(MatchingIndicesFuture[i][2])][int(MatchingIndicesFuture[i][3])]
    MatchingEPWFutureMeanPressure[int(MatchingIndicesFuture[i][0])][int(MatchingIndicesFuture[i][1])] = MeanPressureEPWPast[int(MatchingIndicesFuture[i][2])][int(MatchingIndicesFuture[i][3])]
    MatchingEPWFutureMeanRad[int(MatchingIndicesFuture[i][0])][int(MatchingIndicesFuture[i][1])] = MeanRadEPWPast[int(MatchingIndicesFuture[i][2])][int(MatchingIndicesFuture[i][3])]
    MatchingEPWFutureMeanWindS[int(MatchingIndicesFuture[i][0])][int(MatchingIndicesFuture[i][1])] = MeanWindSEPWPast[int(MatchingIndicesFuture[i][2])][int(MatchingIndicesFuture[i][3])]

    MatchingEPWFutureStdTdryb[int(MatchingIndicesFuture[i][0])][int(MatchingIndicesFuture[i][1])] = StdTdrybEPWPast[int(MatchingIndicesFuture[i][2])][int(MatchingIndicesFuture[i][3])]
    MatchingEPWFutureStdPressure[int(MatchingIndicesFuture[i][0])][int(MatchingIndicesFuture[i][1])] = StdPressureEPWPast[int(MatchingIndicesFuture[i][2])][int(MatchingIndicesFuture[i][3])]
    MatchingEPWFutureStdRad[int(MatchingIndicesFuture[i][0])][int(MatchingIndicesFuture[i][1])] = StdRadEPWPast[int(MatchingIndicesFuture[i][2])][int(MatchingIndicesFuture[i][3])]
    MatchingEPWFutureStdWindS[int(MatchingIndicesFuture[i][0])][int(MatchingIndicesFuture[i][1])] = StdWindSEPWPast[int(MatchingIndicesFuture[i][2])][int(MatchingIndicesFuture[i][3])]

# Use the matching EPW and CM model means and standard deviations to stretch and shift the data for the future period
MatchingEPWFutureTdrybShiftedStretched = numpy.zeros((NYearsFuture, NMonths, NHoursInMonth))
MatchingEPWFuturePressureShiftedStretched = numpy.zeros((NYearsFuture, NMonths, NHoursInMonth))
MatchingEPWFutureRadShiftedStretched = numpy.zeros((NYearsFuture, NMonths, NHoursInMonth))
MatchingEPWFutureWindSShiftedStretched = numpy.zeros((NYearsFuture, NMonths, NHoursInMonth))

for i in range(0, NYearsFuture):
    for j in range(0, NMonths):
        MatchingEPWFutureTdrybShiftedStretched[i][j][:] = MatchingEPWFutureMeanTdryb[i][j] + \
            MeanTdrybCMFuture[i][j] - MatchingEPWFutureMeanTdryb[i][j] + \
            (MatchingEPWFutureTdryb[i][j][:] - MatchingEPWFutureMeanTdryb[i][j]) * \
            StdTdrybCMFuture[i][j] / MatchingEPWFutureStdTdryb[i][j]
        MatchingEPWFuturePressureShiftedStretched[i][j][:] = MatchingEPWFutureMeanPressure[i][j] + \
            MeanPressureCMFuture[i][j] - MatchingEPWFutureMeanPressure[i][j] + \
            (MatchingEPWFuturePressure[i][j][:] - MatchingEPWFutureMeanPressure[i][j]) * \
            StdPressureCMFuture[i][j] / MatchingEPWFutureStdPressure[i][j]
        MatchingEPWFutureRadShiftedStretched[i][j][:] = MatchingEPWFutureMeanRad[i][j] + \
            MeanRadCMFuture[i][j] - MatchingEPWFutureMeanRad[i][j] + \
            (MatchingEPWFutureRad[i][j][:] - MatchingEPWFutureMeanRad[i][j]) * \
            StdRadCMFuture[i][j] / MatchingEPWFutureStdRad[i][j]
        MatchingEPWFutureWindSShiftedStretched[i][j][:] = MatchingEPWFutureMeanWindS[i][j] + \
            MeanWindSCMFuture[i][j] - MatchingEPWFutureMeanWindS[i][j] + \
            (MatchingEPWFutureWindS[i][j][:] - MatchingEPWFutureMeanWindS[i][j]) * \
            StdWindSCMFuture[i][j] / MatchingEPWFutureStdWindS[i][j]

# Check for wind speed not to be negative, otherwise apply a minimum wind speed
for i in range(0, NYearsFuture):
    for j in range(0, NMonths):
        for k in range(0, NHoursInMonth):
            MatchingEPWFutureWindSShiftedStretched[i][j][k] = max(MinimumWindS, MatchingEPWFutureWindSShiftedStretched[i][j][k])

print('\nFinished correcting future data by shifting and stretching.')

# Apply soil temperatures and shift/stretch correction to the future EPW files
# Go to results directory
os.chdir(PathEPWFutureNotCorrected)
# Remember the file names in order in this directory
EPWFutureFiles = os.listdir(PathEPWFutureNotCorrected)

# Adjust EPW files for the future
for i in range(0, NYearsFuture):

    FutureYear = FirstFutureYear + i

    os.chdir(PathEPWFutureNotCorrected)
    EPWSource = EPWFutureFiles[i]
    EPWOutputFileName = EPWMorphedShiftedStretchedFile + numpy.str(FutureYear) + '.epw'
    EPWDestination = PathEPWFutureCorrected + '/' + EPWOutputFileName

    shutil.copyfile(EPWSource, EPWDestination)

    os.chdir(PathEPWFutureCorrected)
    FileDestination = open(EPWDestination, 'r')
    FileDestinationAllLines = FileDestination.readlines()

    for j in range(0, NMonths):

        # Step 1 replace the soil temperature from the matching month

        os.chdir(PathEPWFutureNotCorrected)
        FileDestination = open(EPWDestination, 'w')

        HistoricalYear = int(MatchingIndicesFuture[i * NMonths + j][2])
        HistoricalFileName = EPWHistoricalFiles[HistoricalYear]
        os.chdir(PathEPWHistorical)
        HistoricalFile = open(HistoricalFileName, 'r')

        # Read and discard the first few header lines until soil temperatures are reached
        for k in range(0, SoilTemperatureLineN):
            HistoricalFileLine = HistoricalFile.readline()
            DestinationFileLine = FileDestinationAllLines[k]

        HistoricalSoilTemperatureData = []
        HistoricalSoilTemperatureData.append(list(HistoricalFileLine.split(",")))
        HistoricalSoilTemperatureInput = HistoricalSoilTemperatureData

        DestinationSoilTemperatureData = []
        DestinationSoilTemperatureData.append(list(DestinationFileLine.split(",")))
        DestinationSoilTemperatureInput = DestinationSoilTemperatureData

        for m in range(0, SoilLayersN):
            DestinationSoilTemperatureInput[0][SoilTemperatureStart + m * SoilTemperatureJump + j] = \
                HistoricalSoilTemperatureInput[0][SoilTemperatureStart + m * SoilTemperatureJump + j]

        PrintMe = ""
        for l in range(len(DestinationSoilTemperatureInput[0][:])-1):
            PrintMe += "{}".format(DestinationSoilTemperatureInput[0][l]) + ','
        PrintMe = PrintMe + "{}".format(DestinationSoilTemperatureInput[0][int(len(DestinationSoilTemperatureInput[0][:]))-1])
        NewLine = "{0}".format(PrintMe)

        FileDestinationAllLines[SoilTemperatureLineN-1] = NewLine

        FileDestination.writelines(FileDestinationAllLines)
        FileDestination.close()
        HistoricalFile.close()

        # Step 2: replace corrected weather variables

        os.chdir(PathEPWFutureNotCorrected)
        FileSource = open(EPWSource, 'r')
        FileDestination = open(EPWDestination, 'w')

        # Read and discard the first few header lines plus the previous months
        for k in range(0, NSkipHeader + MonthIndexEPW[j]):
            FileLine = FileSource.readline()

        for k in range(0, MonthIndex[j]):
            FileLine = FileSource.readline()

            ClimateData = []
            ClimateData.append(list(FileLine.split(",")))

            EPWInput = ClimateData

            EPWInput[0][0] = "{0:}".format(FutureYear)
            EPWInput[0][1] = "{0:}".format(j+1)
            # Correct temperature
            EPWInput[0][6] = "{0:.{1}f}".format(
                float(MatchingEPWFutureTdrybShiftedStretched[i][j][k]-CelciusToKelvin), 1)
            # Correct pressure
            EPWInput[0][9] = "{0:.{1}f}".format(
                float(MatchingEPWFuturePressureShiftedStretched[i][j][k]), 1)
            # Correct radiation
            # To replace longwave (12), direct shortwave (14), diffuse shortwave (15) fluxes
            # Calculate a scaling factor based on the total radiation
            RadiationScalingFactor = MatchingEPWFutureRadShiftedStretched[i][j][k] / (
                    float(EPWInput[0][12]) + float(EPWInput[0][14]) + float(EPWInput[0][15]))
            EPWInput[0][13] = "{0:.{1}f}".format(
                float(MatchingEPWFutureRadShiftedStretched[i][j][k]), 1)
            EPWInput[0][12] = "{0:.{1}f}".format(
                float(RadiationScalingFactor * float(EPWInput[0][12])), 1)
            EPWInput[0][14] = "{0:.{1}f}".format(
                float(RadiationScalingFactor * float(EPWInput[0][14])), 1)
            EPWInput[0][15] = "{0:.{1}f}".format(
                float(RadiationScalingFactor * float(EPWInput[0][15])), 1)
            # Correct wind speed
            EPWInput[0][21] = "{0:.{1}f}".format(
                float(MatchingEPWFutureWindSShiftedStretched[i][j][k]), 1)

            PrintMe = ""
            for l in range(len(EPWInput[0][:])-1):
                PrintMe += "{}".format(EPWInput[0][l]) + ','
            PrintMe = PrintMe + "{}".format(EPWInput[0][int(len(EPWInput[0][:]))-1])
            NewLine = "{0}".format(PrintMe)

            FileDestinationAllLines[MonthIndexEPW[j] + k + NSkipHeader] = NewLine

        FileDestination.writelines(FileDestinationAllLines)
        FileDestination.close()
        FileSource.close()

    print('Created corrected EPW file for future year: ', FutureYear)