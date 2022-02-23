#--------------------------------------------------------------------------------------------------------------
# This is a Python3 script used to detrend the annual maximum datas for CONUS
#--------------------------------------------------------------------------------------------------------------

# PIC python3.4.2 has all the following modulars
import numpy as np
import sys, os
import copy
import math
import re
import csv
from shutil import copyfile

####################################################################
UUU = 0;       # 0 if Hist else GCMs
               # change this mannually for different scenario
####################################################################

flag = int(sys.argv[1]); # 0-4 (5 vars)
duration = str(sys.argv[2])  # '24hr', '48hr', '72hr'

if (UUU == 0) :
  scenario = 'Hist'
  print(scenario)
  outdir = '/rcfs/projects/estcp_ngidf/sunn067/detrend_AMF_WY/'+ scenario
  if not os.path.exists(outdir):
    os.makedirs(outdir)
else:
  gcm  = str(sys.argv[3])  # Hist or BNU-ESM
  emiss = str(sys.argv[4]) # emission rcp85 or historical
  scenario = gcm + '/' + emiss
  outdir = './detrend_time_series/'+ gcm
  if not os.path.exists(outdir):
    os.makedirs(outdir)
  outdir = './detrend_time_series/'+ scenario
  if not os.path.exists(outdir):
    os.makedirs(outdir)


outdir = '/rcfs/projects/estcp_ngidf/sunn067/detrend_AMF_WY/%s/%s'  %(scenario, duration)
if not os.path.exists(outdir):
  os.makedirs(outdir)
#--------------------------------------------------------------------------------------------------------------
def npnan(x,y):
	#this function creates the np.nan 2d-array (np.nan should be float)
 	array_2d = np.zeros((x,y), float) 
 	array_2d[:] = np.nan
 	return array_2d

# ---------------------------------------------------------------------
# start the loop, out loop: 3 durations; inner loop: 5 variables, feel free to add more duration or variables
# my expereicen told me the sequential simulation is very fast and donot need running in parallel
# if you want to do it in parallel, there are 2 ways:
# 1) submit multiple scripts concurrently with fixed duratioin or flag
# 2) use cluster rank to assign duration and flag into different rank number


# -------------------------------------------------------------
# identify variable name
if flag == 0:
  variable_name = 'AWR';            trend_var = 'trend_w.txt'
elif flag == 1:
  variable_name = 'max_precip';     trend_var = 'trend_p.txt'
elif flag == 2:
  variable_name = 'max_rainfall';   trend_var = 'trend_r.txt'
elif flag == 3:
  variable_name = 'melt_only';      trend_var = 'trend_melt_o.txt'
else:
  variable_name = 'ROS';            trend_var = 'trend_melt_ros.txt'

outdir = '/rcfs/projects/estcp_ngidf/sunn067/detrend_AMF_WY/%s/%s/%s'  %(scenario, duration, variable_name)
if not os.path.exists(outdir):
 	os.makedirs(outdir)

# -------------------------------------------------------------
# count number of grid (number vary with different GCM)
file_name = './output/trend_summary/%s/%s/%s'  %(scenario, duration, trend_var)  # change the directory for the trend results derived from R script 
lines = [line.rstrip('\n') for line in open(file_name)]	
day = 0
for line in lines:
  item = line.split()
  day += 1
picked_num = day
print('grid number:')
print(picked_num) 


# load the Sen's slope into RAM
sen_slope = npnan(picked_num, 1) 
latt = npnan(picked_num, 1) 
lonn = npnan(picked_num, 1) 

file_name = './output/trend_summary/%s/%s/%s'  %(scenario, duration, trend_var)  # change the directory for the trend results derived from R script 
lines = [line.rstrip('\n') for line in open(file_name)]	
day = 0
for line in lines:
  item = line.split()
  if (item[0] != 'NA'):
    latt[day,0] = float(item[0])  # item[2] indicates the 3rd column in the slope file
    lonn[day,0] = float(item[1])  # item[2] indicates the 3rd column in the slope file
    sen_slope[day,0] = float(item[2])  # item[2] indicates the 3rd column in the slope file
  # the file doesn't exist - no forcing file
  else:
    latt[day,0] = -9999
    lonn[day,0] = -9999
    sen_slope[day,0] = -9999 
  day += 1


# -------------------------------------------------------------
# start detrend
for cell in range(picked_num):
  # detecting no trend, just copying the original file to new directory
  if (sen_slope[cell,0]== 0):   
    from_dir = '/rcfs/projects/estcp_ngidf/sunn067/output_AMF_WY/%s/%s/%s/data_%8.5f_%9.5f.txt'  %(scenario, duration, variable_name, latt[cell,0], lonn[cell,0])
    to_dir   = '/rcfs/projects/estcp_ngidf/sunn067/detrend_AMF_WY/%s/%s/%s/data_%8.5f_%9.5f.txt'  %(scenario, duration, variable_name, latt[cell,0], lonn[cell,0])
    copyfile(from_dir, to_dir)


  # detecting trend
  elif (sen_slope[cell,0]!=-9999):
    # load the original time series into RAM
    file_name = '/rcfs/projects/estcp_ngidf/sunn067/output_AMF_WY/%s/%s/%s/data_%8.5f_%9.5f.txt'  %(scenario, duration, variable_name, latt[cell,0], lonn[cell,0])
    lines = [line.rstrip('\n') for line in open(file_name)]	

    # prepare a np.array to catch the data, assuming the annual max time series has 2 columns: 1) year, 2) annual max values
    var = npnan(len(lines), 2)                # original time series: 0-year; 1-value
    detrend_data = npnan(len(lines), 2)       # detrended time series: 0-year; 1-value

    # load the original data into the var np.array
    day = 0
    for line in lines:
      item = line.split()
      for j in range(len(var[0,:])):
        var[day, j] = float(item[j])
      day += 1

    # copy the 1st column (year)
    detrend_data[:,0] = var[:,0]  


    # detrend the data while keep the mean average
    if (len(lines)%2) == 0:  # even number
      mid_point = int(len(lines)/2 - 1)
      for j in range(0, mid_point):
        detrend_data[j,1] = var[j,1] + sen_slope[cell,0]*(mid_point-j)
      for j in range(mid_point+2, len(lines)):
        detrend_data[j,1] = var[j,1] - sen_slope[cell,0]*(j-mid_point)
      detrend_data[mid_point,1]   = var[mid_point,1]   + sen_slope[cell,0]*0.5
      detrend_data[mid_point+1,1] = var[mid_point+1,1] - sen_slope[cell,0]*0.5
    else:    # odd number
      mid_point = int((len(lines)+1)/2 - 1)
      for j in range(0, mid_point):
        detrend_data[j,1] = var[j,1] + sen_slope[cell,0]*(mid_point-j)
      for j in range(mid_point, len(lines)):
        detrend_data[j,1] = var[j,1] - sen_slope[cell,0]*(j-mid_point)

    # check the mean of detrended values = mean of original values
    try:
      assert(abs(np.mean(detrend_data[:,1])-np.mean(var[:,1]))<2e-1)
    except:
      print('the mean detrended data is not equal to the mean original data')
      print('the duration: %d, the variable name is: %s, the cell number is: %d, the coordinate is: %8.5f, %9.5f') %(duration, variable_name, cell+1, latt[cell,0], lonn[cell,0])
      exit()


    # -------------------------------------------------------------
    # output the file, feel free to change the directory
    outputfile = '/rcfs/projects/estcp_ngidf/sunn067/detrend_AMF_WY/%s/%s/%s/data_%8.5f_%9.5f.txt'  %(scenario, duration, variable_name, latt[cell,0], lonn[cell,0])
    np.savetxt(outputfile, detrend_data, fmt='%d %4.2f')		


  # -------------------------------------------------------------
  # monitoring process
  if  cell==(picked_num-1):
    print('finished processing for scenario: %s, duration: %s and varable: %s' %(scenario, duration, variable_name))
