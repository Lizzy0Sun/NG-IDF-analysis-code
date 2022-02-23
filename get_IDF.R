# -------------------------------------------------------------------------------------------------------
# This is the R script used to estimate the IDF curve for selected return periods
# 2, 5, 10, 25, 50, 100, and 500 years
# the GEV distribution is used here 
# L-moments is used in parameter estimation, do not consider uncertainty at this moment 
# Note: if the data is unable for EVT analaysi (like melt=0 for all years), we put NAN values in IDF
# written by Hongxiang Yan at PNNL for Ning Sun on August 20, 2018
# modified for cluster run by Ning on Aug 31, 2010
# -------------------------------------------------------------------------------------------------------

library(lmom)  # modular written by Hosking, PIC R 3.4.3 has this package

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
  dir.create(sprintf('output/IDF_curves/%s', scenario))
} else {
  gcm = args[2]    # Hist or BNU-ESM/historical
  emiss = args[3]
  scenario <- sprintf('%s/%s', gcm, emiss)
  print(scenario)
  dir.create(sprintf('output/IDF_curves/%s', gcm))
  dir.create(sprintf('output/IDF_curves/%s', scenario))
}
flag = as.integer(args[4])    # 1 - 5
print(flag)


# set the environmental directory, the following relative directory needs to be consistent with it
setwd('/pic/projects/hydronet/CONUS_NGIDF/NGIDF') 

# total num of cells over CONUS (1/16 degree)
num_cell = 207173

# load the coordinates
coor <- data.matrix(read.csv('./dem_livneh', sep = ",", header=FALSE)) # 1-lat; 2-lon; 3-elevation (m)
coor <- coor[,1:2]            # 1st col -> lat; 2nd col -> lon



# -----------------------------------------------------------------------------
# define IDF function, estimate the IDF values
esti_idf <- function(data) {  # data is a n x 1 vector, detrended annual maximum data
  
  # 7 x 1, for each return period
  output <- matrix(data=NaN, nrow=7, ncol=1)  

  # fit the GEV distribution and estimate the parameter using L-moments method
  lmom_bar <- samlmu(data)
  par <- try(pelgev(lmom_bar))
  if (class(par) == 'try-error') { 
    output <- matrix(data=NaN, nrow=7, ncol=1) 
  } else {
      # estimate the quantiles (return periods)
      output[1,1] <- quagev(0.5, para = c(par[[1]], par[[2]], par[[3]]))    # 2-year
      output[2,1] <- quagev(0.8, para = c(par[[1]], par[[2]], par[[3]]))    # 5-year
      output[3,1] <- quagev(0.9, para = c(par[[1]], par[[2]], par[[3]]))    # 10-year
      output[4,1] <- quagev(0.96, para = c(par[[1]], par[[2]], par[[3]]))   # 25-year
      output[5,1] <- quagev(0.98, para = c(par[[1]], par[[2]], par[[3]]))   # 50-year
      output[6,1] <- quagev(0.99, para = c(par[[1]], par[[2]], par[[3]]))   # 100-year
      output[7,1] <- quagev(0.998, para = c(par[[1]], par[[2]], par[[3]]))  # 500-year
  }
  return(output)
}

# -----------------------------------------------------------------------------
# get IDF values for each variable

# loop through each variable and duration
# out loop: variable name
# inner loop: duration
# My experience told me it will take about 5 hours to finish all, so if you want to have results in a hour, you can do it in parallel
# the same like last 2 scripts, there are 2 ways to do in parallel
# 1) fix flag and duration, create 15 scripts and submit 15 jobs concurrently
# 2) use rank num and assign different job into different rank

# range from 1 to 5 variable

  
if (flag == 1) {variable_name = 'AWR'}
if (flag == 2) {variable_name = 'max_precip'}
if (flag == 3) {variable_name = 'max_rainfall'}
if (flag == 4) {variable_name = 'melt_only'}
if (flag == 5) {variable_name = 'ROS'}

# range for 3 duration
for (id in 1:num_cell) {
  
  # load the detrended annual max time series
  file_name <- sprintf('/rcfs/projects/estcp_ngidf/sunn067/detrend_AMF_WY/%s/%s/%s/data_%5.5f_%5.5f.txt', scenario, duration, variable_name, coor[id,1], coor[id,2])
  #file_name <- sprintf('/rcfs/projects/estcp_ngidf/sunn067/detrend_AMF_WY/Hist/24hr/AWR/data_48.84375_-112.65625.txt', scenario, duration, variable_name, coor[id,1], coor[id,2])
  if (file.exists(file_name)) {
  
    # prepare a matrix to store the idf output
    idf_output <- matrix(data=NaN, nrow=7, ncol=1)   # IDF values in mm, the unit is consistent with the annual max values
      
    data <- try(data.matrix(read.table(file_name, header=FALSE)), silent=T)  # 1-col: year; 2-col: annumal max value (mm)
    # if it is an empty file or other reason to cause the try faild, the data returned "try-error"
    if (inherits(data, 'try-error')) {
      idf_output[,1] <- -9999  # assign -9999 to the failed file, 
    } else {
      data <- data[,2]
      # get IDF curves
      output <- esti_idf(data)
      idf_output[,1] <- output[,1]    
    }
    
    
    # monitor the process
    if (id == 1) { cat('start ... \n')}
    if (id == ceiling(num_cell*0.05)) { cat('finished 5% ... \n')}
    if (id == ceiling(num_cell*0.1)) { cat('finished 10% ... \n')}
    if (id == ceiling(num_cell*0.2)) { cat('finished 20% ... \n')}
    if (id == ceiling(num_cell*0.4)) { cat('finished 40% ... \n')}
    if (id == ceiling(num_cell*0.6)) { cat('finished 60% ... \n')}
    if (id == ceiling(num_cell*0.8)) { cat('finished 80% ... \n')}
    if (id == num_cell) { cat(sprintf('finished %s %s %s\n', scenario, variable_name, duration)) }

    # output results into local drive
    dir.create(sprintf('output/IDF_curves/%s/%s', scenario, duration))
    dir.create(sprintf('output/IDF_curves/%s/%s/%s', scenario, duration, variable_name))
    file_name <- sprintf('output/IDF_curves/%s/%s/%s/data_%5.5f_%5.5f.txt', scenario, duration, variable_name, coor[id,1], coor[id,2])
    write.table(idf_output, file_name, row.names=FALSE, col.names=FALSE) 
  }

}














