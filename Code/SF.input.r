
# Initialize Sap Flux Model Inputs

rm(list = ls())

dir <- '/Users/davidbell/StaCC-Model'
setwd(dir)

# Data and Run Name
  DATA   <- 'hardwood'  #'wireless' or 'face'
  RUN    <- 'newcap'    #run name
  TITLE  <- paste(DATA,RUN,sep=".") #prefix to add to output title
	
	# Switched for Optional Features
	  FLAT      <- FALSE    # use flat depth profile?
	  RAND      <- TRUE     # use random effects?
	  BS        <- FALSE    # use canopy status?
	  DPAR      <- TRUE     # estimate data model parameters? 
	  CAP       <- TRUE     # model with capacitance lag
	  moistdef  <- 0.2      #minimum vol moisture for VPD effect estimation
	  daylight  <- FALSE    #if T, only use daylight hours for vpd effect estimation
	  gapfill   <- 3*48     #maximum gap to model
	  Dreg      <- TRUE     #use Gref - lam*ln(D) instead of Gref * (1-lam*ln(D))
	  TC        <-  0.63     #lag for sap flux data
	  SP        <- 'litu'
  
	  intval  <- c(2003.3,2003.8) #select segment to analyze     
 	  year <- trunc(min(intval))
    intval <- intval - year

	  models <- matrix(c(year,SP,TC),ncol=3)
	    colnames(models) <- c('year','SP','TC')
  
	  FLUX.FILE <- "datafiles/js2003.csv"   
	  SAI.FILE  <- paste('datafiles/SA',year,'.csv',sep="")
	  LAI.FILE  <- paste('datafiles/LAI',year,'.csv',sep="")
	  SENS.FILE <- paste('datafiles/HW_JSid.csv',sep="")
	  PAR.FILE  <- paste('datafiles/PAR',year,'.csv',sep="")
	  TEMP.FILE <- paste('datafiles/Ta',year,'.csv',sep="")
	  VPD.FILE  <- paste('datafiles/VPD',year,'.csv',sep="")
	  SM.FILE   <- paste('datafiles/sm',year,'.csv',sep="")
    
#check and select interval
  yrtrunc <- year
  
source('Code/SF.functions.r')
#library(mvtnorm)

options(digits=10)

source('Code/SF.Setup.r')

source('Code/SF.Gibbs.r')


