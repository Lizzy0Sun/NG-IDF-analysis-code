# usage:
#         python3.5 ./construct_idf.py  <stationfile>  <IDF duration>
# Currently called in loop by RunIDF

# purpose:        
# ### build annual maxima time series of following variables:
# 
# ### P:  daily precipitation
# ### R:  rainfall only, when deltaSWE==0
# ### M1: melt only, when deltaSWE<0 and P=0, M1 = -deltaSWE
# ### M2: melt + rain, when deltaSWE<0, M2 = P - deltaSWE
# ### M3: ROS only, when deltaSWE<0 and P > 0, M3 = P - deltaSWE
# ### W: antual amount of water, W = P - deltaSWE all the time


# Import all needed libraries 
import numpy as np
import pandas as pd
import numpy.ma as ma
import math
import sys, os, datetime, time
import re
from pathlib import Path
sys.path.append(os.getcwd())

#--------------------------------------------------------------------------------------------------------------
def npnan(x,y):
    # this function creates the np.nan 2d-array (np.nan should be float)
    array_2d = np.zeros((x,y), float) 
    array_2d[:] = np.nan
    return array_2d 

#--------------------------------------------------------------------------------------------------------------
def finddate(year,month,day,hour,var):
    # this function is given year month day and find the index number in the var array
    # if the year/month/day is not in the var, return the last day of var 
    # input year/month/day: int
    var = var[:,:4]        #only contains the date
    var = var.astype(int)  #float to int
    date = np.array([year,month,day,hour])
    temp = np.where(np.all(var==date, axis=1))  #return a tuple
    try:
        m = int(temp[0])
    except:
        m = len(var[:, 0])-1
    return m

def findMD(a):
    a = a/365*(2*math.pi)
    x = sum(np.cos(a))/len(a)
    y = sum(np.sin(a))/len(a)
    b = np.arctan(y/x)
    
    if (x<0):
        b = math.pi + b
    if (y<0 and x>0):
        b = 2*math.pi + b

    MD = b*365/(2*math.pi)
    return MD

def findSI(a):
    a = a/365*(2*math.pi)
    x = sum(np.cos(a))/len(a)
    y = sum(np.sin(a))/len(a)
    SI = np.sqrt(x*x+y*y)
    return SI

def findDev(a):
    a = a/365*(2*math.pi)
    x = sum(np.cos(a))/len(a)
    y = sum(np.sin(a))/len(a)
    if (x*x+y*y==0):
        dev = 0
    else:
        dev = np.sqrt(-2*np.log(x*x+y*y))*365/(2*math.pi)# deviation
    return dev

def findCV(MD, dev):
    cv = MD/dev
    return cv

maindir = '/pic/projects/hydronet/CONUS_NGIDF/NGIDF/' # where scripts and station files are
outputdir = '/rcfs/projects/estcp_ngidf/sunn067/' # where output will be stored
inputdir = '/rcfs/projects/estcp_ngidf/sunn067/dhsvm_swe_simulations_open/' # where daily SWE simulations are


duration = int(sys.argv[2])  # IDF duration
GCM = str(sys.argv[3])
emiss = str(sys.argv[4])
snotel_num = int(sys.argv[5]) # number of stations in the partitioned station file, e.g. xaa

#--------------------------------------------------------------------------------------------------------------

UUUUU = 0; # 0 is Hist, and 1 if others
if (UUUUU == 1):
    scenario = GCM+'/'+emiss
    file_name = maindir+'stationFiles/'+GCM+'/'+emiss+'/'+str(duration)+'/'+str(sys.argv[1]) # broken statations, e.g. xaa
    #file_name = maindir+'stationFiles/'+str(sys.argv[1]) # broken statations, e.g. xaa
else:
    scenario = 'Hist'
    file_name = maindir+'stationFiles/Hist/'+str(duration)+'/'+str(sys.argv[1]) # broken statations, e.g. xaa
    #file_name = maindir+'stationFiles/'+str(sys.argv[1]) # broken statations, e.g. xaa
print(scenario)

#--------------------------------------------------------------------------------------------------------------
#make output folder if it doesn't exist
if (UUUUU == 1):
    outdir2 = outputdir + 'output_AMF_CY/'+ GCM + '/'
    if not os.path.exists(outdir2):
        os.makedirs(outdir2)
    outdir2 = outputdir + 'output_AMF_CY/'+ GCM + '/' + emiss + '/'
    if not os.path.exists(outdir2):
        os.makedirs(outdir2)
else: 
    outdir2 = outputdir + 'output_AMF_CY/'+ scenario + '/'
#    if not os.path.exists(outdir2):
#        os.makedirs(outdir2)

dir_path1 = outdir2+str(duration)+'hr/max_precip'   
if not os.path.exists(dir_path1):
    os.makedirs(dir_path1)
dir_path2 = outdir2+str(duration)+'hr/max_rainfall'   
dir_path3 = outdir2+str(duration)+'hr/melt_only' 
dir_path4 = outdir2+str(duration)+'hr/AWR'
dir_path5 = outdir2+str(duration)+'hr/ROS'  

#if not os.path.exists(dir_path2):
#    os.makedirs(dir_path2)
  
#if not os.path.exists(dir_path3):
#    os.makedirs(dir_path3)
   
#if not os.path.exists(dir_path4):
#    os.makedirs(dir_path4)            
 
#if not os.path.exists(dir_path5):
#    os.makedirs(dir_path5)


#--------------------------------------------------------------------------------------------------------------
if 'historical' in scenario:
    startyear = 1950        # dhsvm simulation start year
    endyear = 2005          # dhsvm simulation end year
    L1 = 160000
elif 'Hist' in scenario:
    startyear = 1950        # dhsvm simulation start year
    endyear = 2013          # dhsvm simulation end year  
    L1 = 187000
else:
    startyear = 2044        # dhsvm simulation start year
    endyear = 2099          # dhsvm simulation end year
    L1 = 160000    
print(startyear, endyear)

modelStep =3
thresh = 0.2            # ROS threshold

#--------------------------------------------------------------------------------------------------------------
# make time stamps

Tstep = int(duration/modelStep) # reduced array from time aggregation

d1 = datetime.datetime(startyear, 1, 1, 0)  # start date
d2 = datetime.datetime(endyear, 12, 31, 21)  # end date

delta = d2 - d1         # timedelta

dates = []
totalStep = delta.total_seconds()/10800 + 1     # dhsvm output steps
for i in range(int(totalStep-Tstep+1)):
    dates.append(d1 + datetime.timedelta(hours=modelStep*i))
    
    
# read conus grid coords

summary_data = npnan(snotel_num,2)
lines = [line.rstrip('\n') for line in open(file_name)]

count = 0
for line in lines:
    item = line.split(',')
    for i in range(2):
        summary_data[count,i] = float(item[i])
    count += 1



for i in range(snotel_num):    
    
    #(1) load the output data
    dir2 = inputdir+scenario+'/data_%.5f_%.5f' % (summary_data[i,0], summary_data[i,1])
    dir_path = inputdir+scenario+'/data_%.5f_%.5f/WAR.txt' % (summary_data[i,0], summary_data[i,1])
    
    my_file = Path(dir_path)
    count55 = 0
    if (my_file.is_file()):  #and os.stat(dir_path).st_size>100
        for line in open(dir_path).readlines(  ): 
            count55 += 1
    else:
        count55  = 0;
        

    if (count55 > L1) :
        print(i, summary_data[i,0], summary_data[i,1])
        # read files line by line into var1 
        lines = [line.rstrip('\n') for line in open(dir_path)]
        var1 = npnan(len(dates), 5)   # AWR, deltaSWE (t-(t-1)), precipitation, snowfall, swe (all in mm 3-hourly)
        count = 0
        for line in lines[len(lines)-len(dates)-1:-1]:
            item = line.split()
            for j in range(len(var1[0,:])):
                var1[count, j] = float(item[j])
            count += 1
        

        # ## average every duration hours (6 or 24 hours) based on desired duration specified in the beginning.
        var = npnan(int(len(dates)), 9)

        for k in range(int(len(dates))):
            var[k,0] = dates[k].year
            var[k,1] = dates[k].month
            var[k,2] = dates[k].day
            var[k,3] = dates[k].hour

        # average by every duration, take mean over the rpw axis
        for j in range(int(len(dates))):
            for k in range(5):
                if (k==4): # swe: take the last element
                    var[j,k+4] = var1[j:j+Tstep,k][-1]
                else: # take the sum for other vars
                    var[j,k+4] = np.sum(var1[j:j+Tstep,k])
                
                
        # find out the start and end date (index) of each water year
        Startdate = npnan(100, 100000) #start date of each water year
        Enddate = npnan(100, 100000) # end date of each water year 

        k1=0
        k2=0
        for year in np.arange(startyear, endyear, 1):
            for j in range(len(var[:, 0])):
                if (int(var[j, 1]) == 1) and (int(var[j, 0]) == year):
                    Startdate[year-startyear,k1] = j
                    k1 = k1 + 1
                if (int(var[j, 1]) == 12) and (int(var[j, 0]) == year):
                    Enddate[year-startyear,k2] = j
                    k2 = k2 + 1

        startD = []
        endD = []
        for year in np.arange(startyear, endyear, 1):
            startD.append(np.nanmin(Startdate[year-startyear]))
            endD.append(np.nanmax(Enddate[year-startyear]))
        #print(startD, endD)

        ####################################################################################
        #find out the maximum prcp, rainfall, and snowmelt in each year
        ####################################################################################
        max_r  = npnan(100, 3) # rain
        max_m1 = npnan(100, 3) # melt only
        max_m3 = npnan(100, 3) # ROS
        max_w  = npnan(100, 3) # AWR
		max_p  = npnan(100, 3) # precipitation


        count = 0
        for year in np.arange(startyear, endyear, 1):
            
            start_day = int(startD[year-startyear])
            end_day = int(endD[year-endyear])
            
            
            #in the entire water year, the prec and wteq have data for each day
            if (sum(np.isnan(var[start_day:end_day+1, -1])==1)==0) and (sum(np.isnan(var[start_day:end_day+1, 3])==1)==0):
                temp_prcp  = var[start_day:end_day+1, 6]    # precip 
                delta_wteq = var[start_day:end_day+1, 5]    # delta SWE
                temp_w     = var[start_day:end_day+1, 4]    # available amount of water, W = P - deltaSWE               
                temp_r     = var[start_day:end_day+1, 6]-var[start_day:end_day+1, 7]   # R: rainfall only
                temp_r[temp_r<0] = 0;
                swe = var[start_day:end_day+1, 8]    # SWE (mm)

                temp_m1 = npnan(len(temp_prcp), 1)   # melt only, when deltaSWE<0 and P=0, M1 = -deltaSWE  
                temp_m3 = npnan(len(temp_prcp), 1)   # ROS only, when deltaSWE 0, M3 = P - deltaSWE

                #1). find out the max precip 
				max_p[count, 1] = np.nanmax(temp_prcp)    #peak value
				max_p[count, 0] = year+1                  #year
				max_p[count, 2] = int(np.nanargmax(temp_prcp)/8)	
                #--------------------------------------------------------------------------------------------------------
                ### find out the melt only
                for k in range(len(delta_wteq)-1):
                    
                    # M1 (melt only without ROS
                    tot1 = np.abs(delta_wteq[k]) + temp_r[k]
                    if (delta_wteq[k]<0 and temp_r[k]/tot1<thresh and temp_r[k]<(10*duration/24)):
                        temp_m1[k,0] = temp_w[k]


                    # M3 (ROS only) - snowpack of at least 10mm/day SWE, where rain/sum of rainfall and rainfall <=80%
                    if (delta_wteq[k]<0 and swe[k]>=(10*duration/24) and temp_r[k]>=(10*duration/24) and temp_r[k]/tot1<=(1-thresh)): 
                        temp_m3[k,0] = temp_w[k]
                        
                        
                        

                #2). find out the max R1, R2
                max_r[count, 1] = np.nanmax(temp_r)       # peak value
                max_r[count, 2] = int(np.nanargmax(temp_r)/8)    # index of peak value                  
                max_r[count, 0] = year                  # year
                


                #3). find out the max M1, M2, W
                if (sum(np.isnan(temp_m1[:,0]))==len(temp_m1[:,0])):   # if no snow in the entire water year
                    max_m1[count, 1] = 0                               # put 0 value
                    max_m1[count, 2] = -9999
                else: 
                    max_m1[count, 1] = np.nanmax(temp_m1[:,0])         # peak value
                    max_m1[count, 2] = int(np.nanargmax(temp_m1[:,0])/8)    # index of peak value
                max_m1[count, 0] = year                              # year


                if (sum(np.isnan(temp_m3[:,0]))==len(temp_m3[:,0])):   # if no snow in the entire water year
                    max_m3[count, 1] = 0                               # put 0 value
                    max_m3[count, 2] = -9999
                else: 
                    max_m3[count, 1] = np.nanmax(temp_m3[:,0])         # peak value
                    max_m3[count, 2] = int(np.nanargmax(temp_m3[:,0])/8)    # index of peak value
                max_m3[count, 0] = year                              # year


                max_w[count, 1] = np.nanmax(temp_w) ###
                max_w[count, 0] = year
                max_w[count, 2] = int(np.nanargmax(temp_w)/8)

                count += 1
                    
        # summarize mean date, seasonal index, and CV
        
        #remove the extra nan value
		max_p  = max_p[:count,:]
        max_r  = max_r[:count,:]
        max_m1 = max_m1[:count,:]
        max_m3 = max_m3[:count,:]
        max_w  = max_w[:count,:]
        
        # max p                        
        outputfile = dir_path1+'/data_%.5f_%.5f.txt' % (summary_data[i,0], summary_data[i,1])
        np.savetxt(outputfile, max_p, fmt='%d %4.2f %d')
		
		
		# max r
        outputfile = dir_path2+'/data_%8.5f_%9.5f.txt' % (summary_data[i,0], summary_data[i,1])
        np.savetxt(outputfile, max_r, fmt='%d %4.2f %d')


        # m1: melt only
        
        outputfile = dir_path3+'/data_%8.5f_%9.5f.txt' % (summary_data[i,0], summary_data[i,1])
        np.savetxt(outputfile, max_m1, fmt='%d %4.2f %d')


        # w
        outputfile = dir_path4+'/data_%8.5f_%9.5f.txt' % (summary_data[i,0], summary_data[i,1])
        np.savetxt(outputfile, max_w, fmt='%d %4.2f %d')


        # m3: ROS value
        outputfile = dir_path5+'/data_%8.5f_%9.5f.txt' % (summary_data[i,0], summary_data[i,1])
        np.savetxt(outputfile, max_m3, fmt='%d %4.2f %d')

#--------------------------------------------------------------------------------------------------------------
print("safe quit ......")
exit()