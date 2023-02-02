# Morph EPW files for the future period
# Output the files
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
PathEPWFutureNotCorrected = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/EPW-Future-Not-Corrected")
PathResults = ("C:/GoogleDrive/U Guelph/Projects/VWFG/Results/")

MatchFutureFile = "Match-CM-EPW-Future.txt"
EPWTemplateFile = "EPW-Template-Toronto.epw"
#EPWTemplateFile = "EPW-Template-Vancouver.epw"
EPWMorphedFile = "EPW-Morphed-Year-"

# Specify the number of past historical, validation, and future years in the dataset
# Specify the number of months in each year
NYearsPast = 20
NYearsValid = 14
NYearsFuture = 80
NMonths = 12
NHoursInDay = 24
NHoursInMonth = 31 * 24
MatchingFileColumns = 6

CelciusToKelvin = 273.15

# Define starting index for J,F,M,A,M,J,J,A,S,O,N,D in the EPW dataset
# Note: remove February 29 from all the leap years in all the datasets
# Define number of hourly observations in EPW file
MonthIndexCM = [0, 31, 59, 89, 120, 150, 181, 212, 242, 273, 303, 334, 365]
MonthIndexEPW = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
MonthIndex = [744, 672, 744, 720, 744, 720, 744, 744, 720, 744, 720, 744]
EPWNPoints = 8760
DailyNPoints = 365

# Number of rows to skip in the EPW file
NSkipHeader = 8

# Calculations
########################################################################################################################

# Loading and analyzing EPW data for the validation period
# Go to results directory
os.chdir(PathResults)

# Go to directory that contains the historic EPW file
os.chdir(PathEPWHistorical)
# Remember the file names in order in this directory
EPWHistoricalFiles = os.listdir(PathEPWHistorical)

# Loading and analyzing EPW data for the future period
# Go to results directory
os.chdir(PathResults)

# Load matching indices for the future period
MatchingIndicesFuture = numpy.loadtxt(MatchFutureFile)

# Rewrite and create morphed EPW file for the future
for i in range(0, NYearsFuture):

    Source = PathMain + '/' + EPWTemplateFile
    EPWOutputFileName = EPWMorphedFile + numpy.str(i) + '.txt'
    Destination = PathEPWFutureNotCorrected + '/' + EPWOutputFileName

    shutil.copyfile(Source, Destination)

    for j in range(0, NMonths):
        HistoricalYear = int(MatchingIndicesFuture[i * NMonths + j][2])
        YearFile = EPWHistoricalFiles[HistoricalYear]
        StartRowRead = MonthIndexEPW[
            int(MatchingIndicesFuture[i * NMonths + j][3])]      # Row to start reading from the EPW file
        EndRowRead = MonthIndexEPW[
            int(MatchingIndicesFuture[i * NMonths + j][3] + 1)]  # Row to end reading from the EPW file
        StartRowWrite = MonthIndexEPW[j]                         # Row to start writing the EPW file
        EndRowWrite = MonthIndexEPW[j]                           # Row to end writing the EPW file

        os.chdir(PathEPWFutureNotCorrected)
        EPWOutputFile = open(EPWOutputFileName, 'r')
        EPWOutputFileAllLines = EPWOutputFile.readlines()

        os.chdir(PathEPWHistorical)
        File = open(YearFile)

        #Read and discard the first few header lines
        for l in range(0, NSkipHeader):
            FileLine = File.readline()

        # Variable to remember the writing line
        m = NSkipHeader
        for l in range(NSkipHeader, MonthIndexEPW[12] + NSkipHeader + 1):
            FileLine = File.readline()
            if (l >= StartRowRead + NSkipHeader) and (l < EndRowRead + NSkipHeader):
                EPWOutputFileAllLines[StartRowWrite + m] = FileLine
                m = m + 1

        os.chdir(PathEPWFutureNotCorrected)
        EPWOutputFile = open(EPWOutputFileName, 'w')
        EPWOutputFile.writelines(EPWOutputFileAllLines)

        EPWOutputFile.close()

    print('Morphed EPW file for future year: ', i)