#------------------------------------------------------------------------------------------------------------------------
# This is the R script used to do the nonparametric Mann-Kendall test for the annual max data set
# based on the alpha = 5% significance level 
# calculate the trend (Sen's slope) for the annual max time series
# written by Hongxiang Yan at PNNL for Ning Sun, 2018/08/17
#------------------------------------------------------------------------------------------------------------------------

library('trend') # load the M-K test package. Version 1.1.0. Note that the following script does not work for old version
               # It seems PIC does not have trend package. You can try to download to your home directory.
               # The latest "trend" package version is 1.1.1 updated on July 30, 2018, I assumted the following script should work 
               # in the latest version. Let me know if you find any troubles in loading the package.

####################################################################
UUU = 0;       # 0 if Hist else GCMs
               # change this mannually for different scenario
args <- commandArgs(TRUE)
duration = args[1]    # '24hr', '48hr', '72hr' 
print(duration)
####################################################################

if (UUU == 0) {
  scenario = 'Hist'
  print(scenario)
  dir.create(sprintf('trend_summary/%s', scenario))
} else {
  gcm = args[2]    # Hist or BNU-ESM/historical
  emiss = args[3]
  scenario <- sprintf('%s/%s', gcm, emiss)
  print(scenario)
  dir.create(sprintf('trend_summary/%s', gcm))
  dir.create(sprintf('trend_summary/%s', scenario))}


flag = as.integer(args[4])    # 1 - 6
print(flag)


# number of cells, i.e., 207173 for CONUS
num_cell = 207173   

# set the directory and be consistent with the following relative directory such as './'
setwd('/pic/projects/hydronet/CONUS_NGIDF/NGIDF')   
# load the coordinate file for CONUS

coor <- data.matrix(read.csv('./dem_livneh', sep = ",", header=FALSE)) # 1-lat; 2-lon; 3-elevation (m)
coor <- coor[,1:2]

# prepare a matrix to store the trend analysis result
# 1-lat; 2-lon; 3-Sen's slope (unit/yr, unit depending on the annual maximum time series), this value is estimated from Sen's slope. 
# the order in the row number is consistent with the coor file
trend_res <- matrix(data=NaN, nrow=num_cell, ncol=3) 


# -----------------------------------------------------------------------------
# define M-K test function, including the Sen's slope estimation
# this function is assoicated with the "trend" package
# pre-define a R function, which will be used later

MannKendall <- function(data) {   # data is a n x 1 vector
  res <- mk.test(ts(data)) # res is a list of many variables, details can be seen in the R package at: https://cran.r-project.org/web/packages/trend/trend.pdf
  pvalue <- res[[2]] # the p-value for the M-K test
  
  # if all 0 values or the same value (M-K test failed if the data is all 0, or the data series are the same, like 1, 1, 1, 1, 1, ...)
  if (is.nan(pvalue)) {
    slope = 0

  # if pvalue is not NaN (M-K test worked and did not failed)
  } else {  
    if (pvalue > 0.05) {
      slope = 0   # if pvalue > 0.05, there is no significant trend, 0.05 is associated with 5% confidence level
    } else {
      res <- sens.slope(ts(data), conf.level = 0.95)  # if pvalue <= 0.05, there is significant trend and then estimating Sen's slope
      slope = res$estimates
    }
  }
  return(slope)
}

# -----------------------------------------------------------------------------
# M-K test on different valirable 
# range from 1 to 6 variable. This is just for a loop computation. Feel free to revemo or add more variable and change its name.
# let me know if you have any questions on this loop.
# if you think the sequential is too slow, you can do it in parallel. There are two ways:
# 1). create multiple scripts (flag = 1, 2 ,3), and submit the job simultaneously
# 2). use the rank num in the cluster to assign different flag num to different processors (rank num)
# however, my expereicen told me the sequential is fast and do not need parallel

if (flag == 1) {variable_name = 'AWR';          output_name = 'trend_w.txt';}
if (flag == 2) {variable_name = 'max_precip';   output_name = 'trend_p.txt';}
if (flag == 3) {variable_name = 'max_rainfall'; output_name = 'trend_r.txt';}
if (flag == 4) {variable_name = 'melt_only';    output_name = 'trend_melt_o.txt';}
if (flag == 5) {variable_name = 'ROS';          output_name = 'trend_melt_ros.txt';}
if (flag == 6) {variable_name = 'ROS_num';      output_name = 'trend_melt_ros_num.txt';}

  
for (id in 1:num_cell) {
  
  # load the annual max data
  file_name <- sprintf('output/%s/%s/%s/data_%5.5f_%5.5f.txt', scenario, duration, variable_name, coor[id,1], coor[id,2])

  if (file.exists(file_name)) {


    # load the coor
    trend_res[id,1] <- coor[id,1]
    trend_res[id,2] <- coor[id,2]

     # --------------------------------------------------------------------------------------------------------
    # prepare the case when the file is empty. Revision on 09/13/2018
    data <- try(data.matrix(read.table(file_name, header=FALSE)), silent=T)  # 1-col: year; 2-col: annumal max value (mm)
                    
    
    # if it is an empty file or other reason to cause the try faild, the data returned "try-error"
    if (inherits(data, 'try-error')) {
      trend_res[id,3] <- -9999  # assign -9999 to the failed file, 
    } else {
      data <- data[,2]
      # do the M-K
      trend_res[id,3] <- MannKendall(data)      
    }
    # --------------------------------------------------------------------------------------------------------
    
  }
  else {
    trend_res[id,3] <- -9999
  }
    
  # monitor the process
  if (id == 1) { cat('start ... \n')}
  if (id == ceiling(num_cell*0.05)) { cat('finished 5% ... \n')}
  if (id == ceiling(num_cell*0.1)) { cat('finished 10% ... \n')}
  if (id == ceiling(num_cell*0.2)) { cat('finished 20% ... \n')}
  if (id == ceiling(num_cell*0.4)) { cat('finished 40% ... \n')}
  if (id == ceiling(num_cell*0.6)) { cat('finished 60% ... \n')}
  if (id == ceiling(num_cell*0.8)) { cat('finished 80% ... \n')}
}

# output results into local drive
dir.create(sprintf('trend_summary/%s/%s', scenario, duration))
file_name <- sprintf('./trend_summary/%s/%s/%s', scenario, duration, output_name)
write.table(trend_res, file_name, row.names=FALSE, col.names=FALSE) 
    

  
  
  
  

 
  
  
  
  
  
  
















