
#priors
priormatrix <- function(){  #set up prior matrix, values not species-specific

priorlist <- c('prior.A','prior.sdA','loA','hiA','prior.B1','prior.B2',
               'prior.B3','prior.B4',
               'prior.sd1','prior.sd2','prior.sd3','prior.sd4',
               'loB1','loB2','loB3','loB4','hiB1','hiB2','hiB3','hiB4',
               'prior.G','prior.sdG','loG','hiG','prior.ba1','prior.ba2',
               'prior.sda1','prior.sda2',
               'prior.aig',
               'loba1','loba2','hiba1','hiba2',
               'prior.tau','prior.sdtau','lotau','hitau',
			   'prior.beta','prior.sdbeta','lobeta','hibeta',
			   'prior.bstat','lobstat','hibstat','prior.Vbst')
nprior <- length(priorlist)

priormat <- matrix(0,nprior,nspec)
colnames(priormat) <- species
rownames(priormat) <- priorlist

#relative sensitivity to vapor pressure deficit
  #mean
    priormat['prior.A',]    <- .6            
  #standard deviation
    priormat['prior.sdA',]  <- 25            #old   
    if(!Dreg) priormat['prior.sdA',]  <- .2  #new
  #range of reasonable values
    priormat['loA',]        <- 0.4
    priormat['hiA',]        <- 0.8

#prior standard deviation in reference conductance
  priormat['prior.sdG',]  <- 50

#light submodel parameters
  #nighttime conductance (proportion of daytime conductance)
    #mean
      priormat['prior.B1',]   <- .90  
    #standard deviation
      priormat['prior.sd1',]  <- 1
    #range of reasonable values
      priormat['loB1',]       <- .8
      priormat['hiB1',]       <- .99
    
  #sensitivity to light
    #mean
      priormat['prior.B2',]   <- 250 
    #standard deviation
      priormat['prior.sd2',]  <- 100
    #range of reasonable values
      priormat['loB2',]       <- 100
      priormat['hiB2',]       <- 500   

#soil moisture submodel
  #sensitivity to soil moisture
    #mean
      priormat['prior.B3',]   <- .2   #measured moisture
    #standard deviation
      priormat['prior.sd3',]  <- .5
    #range of reasonable values
      priormat['loB3',]       <- .02
      priormat['hiB3',]       <- 0.4

  #threshold soil moisture
    #mean
      priormat['prior.B4',]   <- .2   #max(.25,min(M))
    #standard deviation
      priormat['prior.sd4',]  <- .1
    #range of reasonable values
      priormat['loB4',]       <- 0.15    
      priormat['hiB4',]       <- 0.25

#data model
  #means are species specific
  #standard deviation
    priormat['prior.sda1',]    <- .0001
  #range of reasonable values
    priormat['loba1',]         <- .01  
    priormat['hiba1',]         <- 1

  #means are species specific
  #standard deviation
    priormat['prior.sda2',]    <- .0001
  #range of reasonable values
    priormat['loba2',]         <- .1
    priormat['hiba2',]         <- 1

  #Canopy status effect; i.e., proportional differences between sap flux of dominant vs. supressed trees
    #mean
      priormat['prior.bstat',]   <- 1
    #standard deviation
      priormat['prior.Vbst',]    <- .5
    #range of reasonable values
      priormat['lobstat',]       <- .2
      priormat['hibstat',]       <- 4

#random effects
  priormat['prior.aig',]     <- rep(1,nspec)        #relativized to outer xylem

#range of reasonable values of canopy conductance
  priormat['loG',]       <- 10
  priormat['hiG',]       <- 250

#stomatal lags
  #mean
    priormat['prior.tau',]   <-   10/60/24/365   #based on Woods and Turner 1971 and Naumberg and Ellsworth 2000
  #standard deviation
    priormat['prior.sdtau',] <-   1/60/24/365  
  #range of reasonable values
    priormat['lotau',]       <-    10/60/24/365  
    priormat['hitau',]       <-  10/60/24/365 

#hydraulic lags
  #mean
    priormat['prior.beta',]   <-  .63
  #standard deviation
    priormat['prior.sdbeta',] <-   .01 
  #range of reasonable values
    priormat['lobeta',]       <-  .63
    priormat['hibeta',]       <-  .63

priormat

}

specpriors <- function(sname){  #species-specific priors
                                #for each species, provide the reference conductance (gmean)
                                # and depth profile priors (bamean)

  if(sname %in% c('acru','ACRU')){
    gmean <- 88 
	bamean <- c(.1375,.4555)  #curved acru
              #c(1,.4)      #flat acru
  }
  if(sname %in% c('cato','CATO')){
    gmean <- 44
    bamean <- c(.078,.396)
  }
  if(sname %in% c('LIST','list')){
    gmean <- 85
    bamean <- c(.222,.343)
  }
  if(sname %in% c('LITU','litu')){
    gmean <- 97
    bamean <- c(.055,.568)
  }
  if(sname == 'pita'){
    gmean <- 110
    bamean <- c(1/13.5,0.196) #based on results of Ford et al. 2004
  }
  if(sname %in% c('QUAL','qual')){
    gmean <- 18
    bamean <- c(1,.4)
  }
  if(sname %in% c('QUST','qust')){
    gmean <- 18
    bamean <- c(1,.4)
  }
  if(sname %in% c('QUMI','qumi')){
    gmean <- 18
    bamean <- c(1,.4)
  }
  if(sname %in% c('QUPH','quph')){
    gmean <- 126
    bamean <- c(1,.4)
  }
  if(sname %in% c('QURU','quru')){
    gmean <- 126
    bamean <- c(1,.4)
  }
  
  list(g = gmean, ba = bamean)
}

