
# Clean Workspace
  rm(list = ls())

# Initialize Sap Flux Model Inputs

# Set working directory.
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
	  RELM      <- TRUE     # use relative extractable soil moisture?
	  MULTCAN   <- FALSE    # model as multiple canopies?	
	  CAP       <- TRUE     # model with capacitance lag
	  moistdef  <- 0.2      #minimum vol moisture for VPD effect estimation
	  daylight  <- FALSE    #if T, only use daylight hours for vpd effect estimation
	  gapfill   <- 3*48     #maximum gap to model
	  Dreg      <- TRUE     #use Gref - lam*ln(D) instead of Gref * (1-lam*ln(D))
	  TC        <-  0.63     #lag for sap flux data
	  NC        <- TRUE 
	  SP        <- 'litu'
  
	  intval  <- c(2003.3,2003.8) #select segment to analyze     
 	  year <- trunc(min(intval))
    intval <- intval - year

	  models <- matrix(c(year,SP,TC),ncol=3)
	    colnames(models) <- c('year','SP','TC')
  
	  #Required Data File Locations
	  FLUX.FILE <- paste(dir,'/datafiles/js',year,'.csv',sep="")   
	  SAI.FILE  <- paste(dir,'/datafiles/SA',year,'.csv',sep="")
	  LAI.FILE  <- paste(dir,'/datafiles/LAI',year,'.csv',sep="")
	  SENS.FILE <- paste(dir,'/datafiles/HW_JSid.csv',sep="")
	  PAR.FILE  <- paste(dir,'/datafiles/PAR',year,'.csv',sep="")
	  TEMP.FILE <- paste(dir,'/datafiles/Ta',year,'.csv',sep="")
	  VPD.FILE  <- paste(dir,'/datafiles/VPD',year,'.csv',sep="")
	  SM.FILE   <- paste(dir,'/datafiles/sm',year,'.csv',sep="")
  
#initialize treatments
 
#choose site designation, if not all sites in covariate files are used
  if(!MULTCAN) can.choose <- 'HW'
  
#check and select interval
  yrtrunc <- year
  
source(paste(dir,"Code","SF.functions.r",sep="/"))
#library(mvtnorm)

options(digits=10)

source(paste(dir,"Code","SF.Setup.r",sep="/"))

#boundary line analysis 
  ##need to bin these
 # minQ <- 500
 # minM <- .2
 # bla(minQ, minM)
  
#cov.cor <- cov.plot()

source('SF.Gibbs.r')

}
