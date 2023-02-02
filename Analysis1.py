# Extract data from NetCDF files produced by a CM
# Atmospheric Innovations Research (AIR) Laboratory, University of Guelph
# Developed by Amir A. Aliabadi and Rachel M. McLeod
# Last updated: 2022-03-09

import numpy
import pandas
import netCDF4

# Constants of simulation
########################################################################################################################

# Define latitude and longitude for data extraction
# Note: climate models record a positive longitude, so we add 360 to a negative longitude
lat_rural = 43.649889
lon_rural = 360 - 80.121909

# Define number of years in the dataset
NYears = 5

# Define directories and file names
Directory = "CanRCM4/"
'''
SubDirectory = "NAM-22_ECMWF-ERAINT/"
SubDirectory = "NAM-22_rcp45/"
SubDirectory = "NAM-22_rcp85/"
SubDirectory = "NAM-22_Historical/"
'''

SubDirectory = "NAM-22_rcp45/"

'''
FileNametas = 'tas_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_20010101-20051231.nc'
FileNamesfc = 'sfcWind_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_20010101-20051231.nc'
FileNameps = 'ps_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_20010101-20051231.nc'
FileNamerlds = 'rlds_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_20010101-20051231.nc'
FileNamersds = 'rsds_NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_day_20010101-20051231.nc'
'''

#'''
FileNametas = 'tas_NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_day_20210101-20251231.nc'
FileNamesfc = 'sfcWind_NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_day_20210101-20251231.nc'
FileNameps = 'ps_NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_day_20210101-20251231.nc'
FileNamerlds = 'rlds_NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_day_20210101-20251231.nc'
FileNamersds = 'rsds_NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_day_20210101-20251231.nc'
#'''

'''
FileNametas = 'tas_NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_day_20060101-20091231.nc'
FileNamesfc = 'sfcWind_NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_day_20060101-20091231.nc'
FileNameps = 'ps_NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_day_20060101-20091231.nc'
FileNamerlds = 'rlds_NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_day_20060101-20091231.nc'
FileNamersds = 'rsds_NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_day_20060101-20091231.nc'
'''

tas5Yearnc = Directory + SubDirectory + FileNametas
sfc5Yearnc = Directory + SubDirectory + FileNamesfc
ps5Yearnc = Directory + SubDirectory + FileNameps
rlds5Yearnc = Directory + SubDirectory + FileNamerlds
rsds5Yearnc = Directory + SubDirectory + FileNamersds

'''
OutputFileName0 = Directory + 'NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_2006'
OutputFileName1 = Directory + 'NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_2007'
OutputFileName2 = Directory + 'NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_2008'
OutputFileName3 = Directory + 'NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_2009'
OutputFileName4 = Directory + 'NAM-22_ECMWF-ERAINT_evaluation_r1i1p1_CCCma-CanRCM4_r2_2005'
'''

'''
OutputFileName0 = Directory + 'NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_2001'
OutputFileName1 = Directory + 'NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_2002'
OutputFileName2 = Directory + 'NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_2003'
OutputFileName3 = Directory + 'NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_2004'
OutputFileName4 = Directory + 'NAM-22_CCCma-CanESM2_historical_r1i1p1_CCCma-CanRCM4_r2_2005'
'''

#'''
OutputFileName0 = Directory + 'NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_2021'
OutputFileName1 = Directory + 'NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_2022'
OutputFileName2 = Directory + 'NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_2023'
OutputFileName3 = Directory + 'NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_2024'
OutputFileName4 = Directory + 'NAM-22_CCCma-CanESM2_rcp45_r1i1p1_CCCma-CanRCM4_r2_2025'
#'''

# Calculations and writing data to file
########################################################################################################################

# Load the dataset
tas5Year = netCDF4.Dataset(tas5Yearnc)
sfc5Year = netCDF4.Dataset(sfc5Yearnc)
ps5Year = netCDF4.Dataset(ps5Yearnc)
rlds5Year = netCDF4.Dataset(rlds5Yearnc)
rsds5Year = netCDF4.Dataset(rsds5Yearnc)

# Shows the content of the dataset
print(tas5Year.variables.keys())
print(sfc5Year.variables.keys())
print(ps5Year.variables.keys())
print(rlds5Year.variables.keys())
print(rsds5Year.variables.keys())

# Access variables inside NetCDF file
timetas5Year = tas5Year.variables['time'][:]
timesfc5Year = sfc5Year.variables['time'][:]
timeps5Year = ps5Year.variables['time'][:]
timerlds5Year = rlds5Year.variables['time'][:]
timersds5Year = rsds5Year.variables['time'][:]

lattas5Year = tas5Year.variables['lat'][:]
latsfc5Year = sfc5Year.variables['lat'][:]
latps5Year = ps5Year.variables['lat'][:]
latrlds5Year = rlds5Year.variables['lat'][:]
latrsds5Year = rsds5Year.variables['lat'][:]

lontas5Year = tas5Year.variables['lon'][:]
lonsfc5Year = sfc5Year.variables['lon'][:]
lonps5Year = ps5Year.variables['lon'][:]
lonrlds5Year = rlds5Year.variables['lon'][:]
lonrsds5Year = rsds5Year.variables['lon'][:]

tas5Year = tas5Year.variables['tas'][:]
sfc5Year = sfc5Year.variables['sfcWind'][:]
ps5Year = ps5Year.variables['ps'][:]
rlds5Year = rlds5Year.variables['rlds'][:]
rsds5Year = rsds5Year.variables['rsds'][:]

print('Dimension of time:', timetas5Year.shape)
print('Dimension of latitude and longitude:', lattas5Year.shape)
print('Dimension of weather variables:', tas5Year.shape)

LatDiff = 1e6
LonDiff = 1e6
for i in range(0, lattas5Year.shape[0]):
    for j in range(0, lattas5Year.shape[1]):
        if numpy.abs(lattas5Year[i][j] - lat_rural) < LatDiff:
            LatDiff = numpy.abs(lattas5Year[i][j] - lat_rural)
            iIndex  = i
        if numpy.abs(lontas5Year[i][j] - lon_rural) < LonDiff:
            LonDiff = numpy.abs(lontas5Year[i][j] - lon_rural)
            jIndex  = j

print('iIndex, jIndex, Lat Match, lon Match [deg] = ', iIndex, jIndex, lattas5Year[iIndex][jIndex], lontas5Year[iIndex][jIndex])

NTimeStamps = len(timetas5Year)
NTimeStampsInYear = int(NTimeStamps/NYears)
print(NTimeStamps)

OutputFile0 = open(OutputFileName0, "w")
OutputFile1 = open(OutputFileName1, "w")
OutputFile2 = open(OutputFileName2, "w")
OutputFile3 = open(OutputFileName3, "w")
OutputFile4 = open(OutputFileName4, "w")

OutputFile0.write("# 0: Day of Year [-] \t 1: lat [deg] \t 2: lon [deg] \t 3: tas [K] \t 4: sfcWind [m s-1] \t 5: ps [Pa] \t 6: Total Downwelling Rad Flux [W m-2] \n")
OutputFile1.write("# 0: Day of Year [-] \t 1: lat [deg] \t 2: lon [deg] \t 3: tas [K] \t 4: sfcWind [m s-1] \t 5: ps [Pa] \t 6: Total Downwelling Rad Flux [W m-2] \n")
OutputFile2.write("# 0: Day of Year [-] \t 1: lat [deg] \t 2: lon [deg] \t 3: tas [K] \t 4: sfcWind [m s-1] \t 5: ps [Pa] \t 6: Total Downwelling Rad Flux [W m-2] \n")
OutputFile3.write("# 0: Day of Year [-] \t 1: lat [deg] \t 2: lon [deg] \t 3: tas [K] \t 4: sfcWind [m s-1] \t 5: ps [Pa] \t 6: Total Downwelling Rad Flux [W m-2] \n")
OutputFile4.write("# 0: Day of Year [-] \t 1: lat [deg] \t 2: lon [deg] \t 3: tas [K] \t 4: sfcWind [m s-1] \t 5: ps [Pa] \t 6: Total Downwelling Rad Flux [W m-2] \n")


for k in range(0, NTimeStampsInYear):
    OutputFile0.write("%i \t %f \t %f \t %f \t %f \t %f \t %f \n"
                          % (k, lattas5Year[iIndex][jIndex], lontas5Year[iIndex][jIndex], tas5Year[k][iIndex][jIndex],
                             sfc5Year[k][iIndex][jIndex], ps5Year[k][iIndex][jIndex],
                          rlds5Year[k][iIndex][jIndex] + rsds5Year[k][iIndex][jIndex]))
    OutputFile1.write("%i \t %f \t %f \t %f \t %f \t %f \t %f \n"
                          % (k, lattas5Year[iIndex][jIndex], lontas5Year[iIndex][jIndex], tas5Year[k + NTimeStampsInYear][iIndex][jIndex],
                             sfc5Year[k + NTimeStampsInYear][iIndex][jIndex], ps5Year[k + NTimeStampsInYear][iIndex][jIndex],
                          rlds5Year[k + NTimeStampsInYear][iIndex][jIndex] + rsds5Year[k + NTimeStampsInYear][iIndex][jIndex]))
    OutputFile2.write("%i \t %f \t %f \t %f \t %f \t %f \t %f \n"
                          % (k, lattas5Year[iIndex][jIndex], lontas5Year[iIndex][jIndex], tas5Year[k + 2 * NTimeStampsInYear][iIndex][jIndex],
                             sfc5Year[k + 2 * NTimeStampsInYear][iIndex][jIndex], ps5Year[k + 2 * NTimeStampsInYear][iIndex][jIndex],
                          rlds5Year[k + 2 * NTimeStampsInYear][iIndex][jIndex] + rsds5Year[k + 2 * NTimeStampsInYear][iIndex][jIndex]))
    OutputFile3.write("%i \t %f \t %f \t %f \t %f \t %f \t %f \n"
                          % (k, lattas5Year[iIndex][jIndex], lontas5Year[iIndex][jIndex], tas5Year[k + 3 * NTimeStampsInYear][iIndex][jIndex],
                             sfc5Year[k + 3 * NTimeStampsInYear][iIndex][jIndex], ps5Year[k + 3 * NTimeStampsInYear][iIndex][jIndex],
                          rlds5Year[k + 3 * NTimeStampsInYear][iIndex][jIndex] + rsds5Year[k + 3 * NTimeStampsInYear][iIndex][jIndex]))
    OutputFile4.write("%i \t %f \t %f \t %f \t %f \t %f \t %f \n"
                          % (k, lattas5Year[iIndex][jIndex], lontas5Year[iIndex][jIndex], tas5Year[k + 4 * NTimeStampsInYear][iIndex][jIndex],
                             sfc5Year[k + 4 * NTimeStampsInYear][iIndex][jIndex], ps5Year[k + 4 * NTimeStampsInYear][iIndex][jIndex],
                          rlds5Year[k + 4 * NTimeStampsInYear][iIndex][jIndex] + rsds5Year[k + 4 * NTimeStampsInYear][iIndex][jIndex]))

OutputFile0.close()
OutputFile1.close()
OutputFile2.close()
OutputFile3.close()
OutputFile4.close()



