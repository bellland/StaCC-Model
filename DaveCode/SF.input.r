
# Initialize Sap Flux Model Inputs

rm(list = ls())

setwd('/Volumes/dbell9$/Sap Flux')

# Data and Run Name
  DATA   <- 'hardwood'  #'wireless' or 'face'
  RUN    <- 'newcap'    #run name
  TITLE  <- paste(DATA,RUN,sep=".") #prefix to add to output title

#
  if(DATA == 'hardwood'){ #begin hardwood

	#for cluster run automation REMOTE==TRUE
    	spCluster <- c('list','qumi') #c('cato','list','litu','qual','qumi','quph') #species to run
    	tcCluster <- c(0.63) #lag parameters
	
	# Switched for Optional Features
	  FLAT      <- FALSE    # use flat depth profile?
	  SECTION   <- FALSE    # T = use three sapwood sections F = use continuous depth profile
	  RAND      <- TRUE     # use random effects?
	  BS        <- FALSE    # use canopy status?
	  DPAR      <- TRUE     # estimate data model parameters? 
	  RELM      <- TRUE     # use relative extractable soil moisture?
	  MULTCAN   <- FALSE    # model as multiple canopies?	
	  CAP       <- TRUE     # model with capacitance lag
	  moistdef  <- 0.2      #minimum vol moisture for VPD effect estimation
	  daylight  <- FALSE    #if T, only use daylight hours for vpd effect estimation
	  REMOTE    <- FALSE    # if on cluster and you want to run each year seperately, set to TRUE         
	  REPEAT    <- FALSE    # if not REMOTE, but one still wishes to do runs across years, species, and lags
	  gapfill   <- 3*48     #maximum gap to model
	  Dreg      <- TRUE     #use Gref - lam*ln(D) instead of Gref * (1-lam*ln(D))
	  MultYear  <- FALSE    #Run years simultaneously
	  TC        <-  0.63     #lag for sap flux data
	  NC        <- TRUE 
	  SP        <- 'litu'
  
	  intval  <- c(2003.3,2003.8) #select segment to analyze     
 	  year <- trunc(min(intval))

	  if(MultYear) year   <- c(2002,2003,2004,2005)

	  models <- matrix(c(year,SP,TC),ncol=3)
	    colnames(models) <- c('year','SP','TC')
  
	  if(REPEAT) {
	  	models <- expand.grid(year=c(2002,2003,2004,2005),SP=spCluster,TC=tcCluster)

	  	models <- models[models[,'SP'] != 'qumi' | models[,'year']>=2004,]
	  
	  }
  
	  runs <- 1 #(1:nrow(models))

	  for(j in runs){
	  
	  year <- as.integer(models[j,'year'])
	  SP   <- as.character(models[j,'SP'])
	  TC   <- as.numeric(models[j,'TC'])
  

		if(REMOTE) {
			tmp <- expand.grid(c(2002,2003,2004,2005),spCluster,tcCluster)
			year <- tmp[as.integer(Sys.getenv('SGE_TASK_ID')),1]
	    	SP <- as.character(tmp[as.integer(Sys.getenv('SGE_TASK_ID')),2])
	    	TC <- tmp[as.integer(Sys.getenv('SGE_TASK_ID')),3]
		}
		
		
		  intval  <- c(.3,.8) #select segment to analyze     

	  #Required Data File Locations
	if(!MultYear){
	  FLUX.FILE <- paste('datafiles/js',year,'.csv',sep="")   
	  SAI.FILE  <- paste('datafiles/SA',year,'.csv',sep="")
	  LAI.FILE  <- paste('datafiles/LAI',year,'.csv',sep="")
	  SENS.FILE <- 'datafiles/HW_JSid.csv'
	  PAR.FILE  <- paste('datafiles/PAR',year,'.csv',sep="")
	  TEMP.FILE <- paste('datafiles/Ta',year,'.csv',sep="")
	  VPD.FILE  <- paste('datafiles/VPD',year,'.csv',sep="")
	  SM.FILE   <- paste('datafiles/sm',year,'.csv',sep="")
  

	#Optional Data File Locations
	  if(SECTION) DIST.FILE <- 'datafiles/ASdist.csv'
	}

	if(MultYear){
	  FLUX.FILE <- paste('datafiles/js',year,'.csv',sep="")   
	  SAI.FILE  <- paste('datafiles/SA',year,'.csv',sep="")
	  LAI.FILE  <- paste('datafiles/LAI',year,'.csv',sep="")
	  SENS.FILE <- 'datafiles/HW_JSid.csv'
	  PAR.FILE  <- paste('datafiles/PAR',year,'.csv',sep="")
	  TEMP.FILE <- paste('datafiles/Ta',year,'.csv',sep="")
	  VPD.FILE  <- paste('datafiles/VPD',year,'.csv',sep="")
	  SM.FILE   <- paste('datafiles/sm',year,'.csv',sep="")
	  

	#Optional Data File Locations
	  if(SECTION) DIST.FILE <- 'datafiles/ASdist.csv'
	}
  
	}

	} #end hardwood
  
#initialize treatments
  ntreat <- 0
  nfact  <- 0
  if(DATA == 'face'){
    treats <- c('Cb','Ft')
    
	Cb<-c('A','E')
    Ft<-c('C','F')
	ii    <- 1
	year <- 2006
	  #check and select interval
  intval  <- c(2006,2006.5) #select segment to analyze    
  year <- trunc(min(intval))
  intval <- intval - year
 
  treatmat <- expand.grid(get(treats[1]),get(treats[2]))
	  treatmat <- matrix(as.character(unlist(treatmat)),nrow=dim(treatmat)[1])
      colnames(treatmat) <- treats
	  
    treatmat <- cbind(treatmat,Treat=paste(treatmat[,1],treatmat[,2],sep=""))
	nfact  <- length(treats)
	ntreat <- dim(treatmat)[1]

	Treat<-treatmat[ii,'Treat']
    TreatPeriod<-paste(treatmat[ii,'Treat'],year,sep='')
    if(year<2001|year>2004|treatmat[ii,'Ft']=='C'){
    if(treatmat[ii,'Cb']=='A') rings<-cbind("A",c(1,5,6,8)) else rings<-cbind('E',c(2,3,4,7))}else{
	  if(treatmat[ii,'Cb']=='A') rings<-cbind("A",c(8)) else rings<-cbind("E",c(7))
	  }

    Fert<-treatmat[ii,'Ft']
    Carb<-treatmat[ii,'Cb']
    SENS.FILE <-paste('datafiles/',Carb,Fert,year,'cleanJSid.csv',sep='')
    FLUX.FILE <- paste('datafiles/',Carb,Fert,year,'cleanJS.csv',sep='')
    PAR.FILE  <-paste('datafiles/','PAR',year,'rescaled',Carb,'.csv',sep='')
    SM.FILE   <-paste('datafiles/','M',year,Carb,'.csv',sep='')
    VPD.FILE  <-paste('datafiles/','AveUpD',year,Carb,'.csv',sep='')
    TEMP.FILE <- paste('datafiles/','AirTempUp',Carb,'.csv',sep='')
    SAI.FILE  <- 'datafiles/SAInew.csv'
    LAI.FILE  <- 'datafiles/DailyLAInew.csv'
    if(SECTION) DIST.FILE <- 'datafiles/ASdist.csv'

  }


#choose site designation, if not all sites in covariate files are used
  if(!MULTCAN) can.choose <- 'HW'
  
#check and select interval
  yrtrunc <- year
  
source('SF.functions.r')
#library(mvtnorm)

options(digits=10)

source('SF.Setup.r')

#boundary line analysis 
  ##need to bin these
 # minQ <- 500
 # minM <- .2
 # bla(minQ, minM)
  
#cov.cor <- cov.plot()

source('SF.Gibbs.r')

}
