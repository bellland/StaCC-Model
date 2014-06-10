# Import and Organize Data for Sap Flux Model

# Import Data
tmp <- inData()
 Q          <- tmp$Q                #light data
 M          <- tmp$M                #soil moisture data
 D          <- tmp$D                #vpd data
 Temp       <- tmp$Temp             #temperature data
 ntree      <- tmp$ntree            #number of trees
 nprobe     <- tmp$nprobe           #number of probes
 depth      <- tmp$depth            #probe depth
 Status     <- tmp$Status           #tree canopy status
 Jdata      <- tmp$Jdata            #sap flux data
 starts     <- tmp$starts           #start point for each time series
 stops      <- tmp$stops            #end point for each time series
 gaps       <- tmp$gaps             #
 timeall    <- tmp$timeall                      #time
 dt         <- diff(timeall)                    #time intervals between time steps
 nt         <- tmp$nt                           #number of timesteps
 LAIspecies <- tmp$LAIspecies                   #leaf area for each species
 LAItree    <- tmp$LAItree                      #leaf area index for each probe
 specindex  <- tmp$specindex                    #species ID for each probe
 SAPtree    <- tmp$SAPtree                      #sapwood area for each probe
 SAPspecies <- tmp$SAPspecies                   #sapwood area for each probe
 if(SECTION)SAratio    <- as.numeric(tmp$SAratio[,4:dim(tmp$SAratio)[2]])  #distribution of sapwood when SECTION == T
 probe      <- tmp$probe                        #probe names
 species    <- tmp$species                      #list of species
 nspec      <- tmp$nspec                        #number of species
 specvars   <- tmp$specvars                     #table of leaf and sapwood areas for each species at each site
 rad        <- tmp$rad
 #scale of parameter estimation
   SITE       <- tmp$SITE                         #site IDs
   nsite      <- length(SITE)                    #number of sites
   siteindex  <- tmp$siteindex                    #site ID for each probe
 #scale of covariates
   plots      <- tmp$plots                        #plot IDs
   nplot      <- length(plots)                   #number ofplots 
   plotindex  <- tmp$plotindex                    #plot ID for each probe
  
#restrict analysis to range of dates provided in intval 
 tmp <- getSection()     
    starts <- tmp$starts
     stops <- tmp$stops
     Jdata <- tmp$Jdata
      Temp <- tmp$Temp
         Q <- tmp$Q
         M <- tmp$M
         D <- tmp$D
   LAItree <- tmp$LAItree
LAIspecies <- tmp$LAIspecies
   SAPtree <- tmp$SAPtree
SAPspecies <- tmp$SAPspecies
   timeall <- tmp$timeall
        dt <- tmp$dt
        nt <- length(timeall)


 D[D<=0] <- 0.001
 deficit_0 <- 0
 deficit <- matrix(deficit_0,nt,nsite)
 wdD     <-  which(D < .6,arr.ind=T)

#remove columns with no data
    keep   <- which(colSums(Jdata,na.rm=T)>0)
   Jdata   <- Jdata[,keep]
   rad     <- rad[keep,]
 notmiss   <- which(is.finite(Jdata),arr.ind=T)
   ntree   <- length(keep)
  nprobe   <- ntree
   depth   <- depth[keep]
  Status   <- 1 - Status[keep]     #changes status to 1=suppressed, 0=dominant
 LAItree   <- matrix(LAItree[,keep],nt)
 SAPtree   <- matrix(SAPtree[,keep],nt)
 specindex <- species[specindex[keep]]
 species   <- sort(unique(specindex))
 specindex <- match(specindex,species)
 siteindex <- siteindex[keep]
 plotindex <- plotindex[keep]
 #numeric place holders for sites and species
 site.spec  <- cbind(siteindex,specindex)
   colnames(site.spec) <- c('SITE','specindex')
 site.all   <- matrix(unique(site.spec),ncol=2)
 site.all <- matrix(site.all[order(paste(site.all[,1],site.all[,2])),],ncol=2)
   colnames(site.all) <- c('SITE','specindex')
 site.sensor <- match(paste(site.spec[,1],site.spec[,2]),
                      paste(site.all[,1],site.all[,2]))
 LAIspecies<- matrix(LAIspecies[,sort(unique(site.sensor))],nt)
 SAPspecies<- matrix(SAPspecies[,sort(unique(site.sensor))],nt)

 probe     <- probe[keep]
 nspec     <- length(species)

 DEPTH      <- sort(unique(depth))
 depthindex <- match(depth,DEPTH)
 
 #priors
 
priormat <- priormatrix()

for(j in 1:nspec){
  ptmp <- specpriors(species[j])
  priormat['prior.G',j]   <- ptmp$g
  priormat['prior.ba1',j] <- ptmp$ba[1]
  priormat['prior.ba2',j] <- ptmp$ba[2]
}

  if (FLAT) {
  priormat['prior.ba1',] <- 1.0
  priormat['prior.ba2',] <- 0.4
  priormat['prior.aig',] <- 1.0
  }
  
if(!DPAR) {
  priormat['loba1',]      <- priormat['prior.ba1',]
  priormat['loba2',]      <- priormat['prior.ba2',]
  priormat['hiba1',]      <- priormat['prior.ba1',]
  priormat['hiba2',]      <- priormat['prior.ba2',]
  }
  
if (!BS) {
  priormat['lobstat',] <- 1
  priormat['hibstat',] <- 1
  }
if (BS) {
  tmp <- tapply(apply(Jdata,2,mean,na.rm=T),paste(specindex,Status),mean)
  tmp <- tmp[(1:nspec)*2]/tmp[(1:nspec)*2-1]
  priormat['prior.bstat',] <- tmp 
  }

#use nplot or nsite for priors?  
#prior for process error
mus <- 5
s1 <- nt*nspec*nsite
s2 <- mus*(s1 - 1)

#prior for measurement error
muv <- 5
v1 <- length(which(is.finite(Jdata)))/10
v2 <- muv*(v1 - 1)

#prior for random effects error
meana <- priormat['prior.aig',] 
mva <- 1
va1 <- max(nprobe,2)  #priors for variance for random effects
va2 <- mva*(va1 - 1)

if(SECTION){
  vb <- .2
  #mvb <- .2
  #vb1 <- nprobe/3
  #vb2 <- mvb*(vb1 - 1)
  #vb <- mvb
}
#initial parameter values

gspec      <- priormat['prior.G',]
aspec      <- priormat['prior.A',]
if(Dreg) 
     aspec <- gspec*aspec

minpropa <- priormat['loA',1]
maxpropa <- priormat['hiA',1]

blitespec  <- matrix(priormat[c('prior.B1','prior.B2'),],nrow=2)
bmoistspec  <- matrix(priormat[c('prior.B3','prior.B4'),],nrow=2)
bag        <- matrix(priormat[c('prior.ba1','prior.ba2'),],nrow=2)

	for(j in 1:nspec)
	  for(i in 1:nsite){
	    specvars[species[j],paste('SAP',SITE[i],sep="-")] <- 
		  sum(pi*rad[specindex==j & siteindex==i,2]^2 - pi*rad[specindex==j & siteindex==i,1]^2)/10000
		specvars[species[j],paste('n',SITE[i],sep="-")] <- dim(rad)[1]
		}

		  
if(SECTION) {
  bag <- matrix(1,length(DEPTH),nspec)
  for(k in 1:nspec)
    for(j in 1:length(DEPTH))
      bag[j,] <- mean(Jdata[,specindex==k & depthindex==j],na.rm=T)/mean(Jdata[,specindex==k],na.rm=T)
  
    prior.ba <- bag
	loba <- bag*0
	hiba <- bag*2

  if(!DPAR){
    prior.ba <- bag
    loba      <- bag
    hiba      <- bag
    }
  }

if(!SECTION) {   #rescale?
    prior.ba <- bag
  loba <- rbind(priormat['loba1',],priormat['loba2',])
  hiba <- rbind(priormat['hiba1',],priormat['hiba2',])
  bag  <- rbind(priormat['prior.ba1',],priormat['prior.ba2',])
  sb1 <- priormat['prior.sda1',]
  sb2 <- priormat['prior.sda2',]
  SPECVARS <- specvars
  tmp <- rescale()
    specvars   <- tmp$sv2
  } 

priorgt.lo <- matrix(rep(0,nt),nt,dim(site.all)[1])
priorgt.hi <- matrix(rep(500,nt),nt,dim(site.all)[1])
priorgt.hi[Q[,1] < 1,] <- .2*500  #nighttime conductance can be 20% of maximum

Gsmat <- gssMat(gspec,aspec,blitespec,bmoistspec)
  Gsmat[Gsmat>priorgt.hi] <- priorgt.hi[Gsmat>priorgt.hi]
  Gsmat[Gsmat<priorgt.lo] <- priorgt.lo[Gsmat<priorgt.lo]
 Gtmat <- Gsmat

tau    <- priormat['prior.tau',]            
beta  <- priormat['prior.beta',]            
if(!CAP) beta <- beta*0 + min(dt)
verror <- runif(1,0,20) 
sigma  <- runif(1,0,20)
werror <- runif(1,0,.1)

bstat <- priormat['prior.bstat',]
if (!BS) bstat <- rep(1,nspec)

aig <- priormat['prior.aig',specindex]  #random intercepts for data model
lor <- aig                          #lowest random effect value
hir <- aig                          #highest random effect value
if(RAND) {
  lor <- 0.25                           #25% of the minimum mean intercept (1)
  hir <- 4.00                           #400% of the maximum mean intercept (1.5)
}
ba  <- bag                              #mean, sd for data model
  
Jtmat <- matrix(NA,nrow=nt,ncol=dim(site.all)[1])
e.qt <- matrix(qt(LAIspecies,SAPspecies,D[,site.all[,1]],Temp[,site.all[,1]]),nt)
for(j in 1:dim(site.all)[1]) {
  if(length(which(site.sensor==j))==1)
    Jtmat[,j] <- Jdata[,site.sensor==j]
  if(length(which(site.sensor==j))>1)
    Jtmat[,j] <- apply(Jdata[,site.sensor==j],1,mean,na.rm=T)
  }

nomean <- which(is.na(Jtmat),arr.ind=T)
  if(length(nomean)>0)
    Jtmat[nomean] <- Gtmat[nomean]*e.qt[nomean]
Jpred<-Jtmat

#initialize the water storage matrix
Wtmat <- Jtmat*0

  if(!year %in% seq(1980,2020,by=4))
     DT <- dt*365*24*60*60
  
  if(year %in% seq(1980,2020,by=4))
     DT <- dt*366*24*60*60


for(t in 2:nt) {
  tmp <- WJcalc(Wtmat[t-1,],Gtmat[t-1,]*e.qt[t-1,],DT[t-1],beta)
    Wtmat[t,] <- tmp$W
    Jtmat[t,] <- tmp$J
  }


if(CAP) gtmp    <- update_gt_CAP()
if(!CAP) gtmp    <- update_gt_NOCAP()
Gtmat   <- gtmp$Gtmat

if(CAP){
  Jtmat   <- gtmp$Jtmat
  Wtmat   <- gtmp$Wtmat
  }
if(!CAP)
  Jtmat <- Gtmat*e.qt
  
Jpred <- matrix(0,nt,nprobe)
pj <- pred_jt(aig,bag,bstat,SECTION)
  Jpred   <- Jtmat[,site.sensor]*matrix(pj,nt,nprobe,byrow=T)
  colnames(Jpred) <- probe
  

pc     <- rbind(c(.005,-.0001),c(-.0001,.005))
if(SECTION) pc <- diag(rep(.1,length(bag)))
pcovba <- kronecker(diag(1,nspec),pc)         #data parameters
Ib    <- pcovba
Ib[which(Ib!=0)] <- 1

if(CAP) pcovK <- diag(1,nspec)

pcovbs <- rep(.0002,nspec)

nd <- nt - dim(wdD)[1]
pcova  <- rbind(c(1,-.05),c(-.05,.01))
if(Dreg) pcova <- rbind(c(1,-.05),c(-.05,1))
pcovga <- kronecker(pcova,diag(1,nspec))
Ig <- pcovga
Ig[Ig!=0] <- 1

pc    <- rbind(c(.01,0),c(0,10))
pcovQ <- kronecker(diag(1,nspec),pc)
Ik    <- pcovQ
Ik[which(Ik!=0)] <- 1

pc    <- rbind(c(.01,-.00001),c(-.00001,.01))
pcovM <- kronecker(diag(1,nspec),pc)
#pcovM  <- rep(.08,nspec)
va <- 1

  pcovb  <- diag(c(1,.1,.01,.01))
 # pcovb  <- solve(crossprod(cbind(rep(1,nt),Q,M,Q*M)))
 

ng    <- 10000
  gburnin <- 3000 
kg    <- seq(40,ng,by=40)
saveg <- ng
if(ng > 2000) saveg <- seq(200,ng,by=200)

dgibbs  <- matrix(NA,ng,length(bag) + nspec)
colnames(dgibbs) <- c(as.vector(outer(paste('b',1:dim(bag)[1],sep=''),
                      species,paste,sep='-')),paste('bstat',species,sep='-'))
agibbs  <- matrix(0,ng,nprobe)
lgibbs  <- matrix(0,ng,length(blitespec))
colnames(lgibbs) <- as.vector(outer(c('l1','l2'),species,paste,sep='-'))
mgibbs  <- matrix(0,ng,length(bmoistspec))
colnames(mgibbs) <- as.vector(outer(c('m1','m2'),species,paste,sep='-'))
ggibbs <- matrix(0,ng,nspec)
colnames(ggibbs) <- species
agibbs <- ggibbs

vgibbs  <- matrix(NA,ng,3+2*nspec)
colnames(vgibbs) <- c('sigma','verror','va',paste('tau',species,sep="-"),paste('beta',species,sep="-"))
rgibbs  <- matrix(NA,ng,nprobe)
colnames(rgibbs) <-probe 
Ggibbs  <- matrix(0,nt,dim(site.all)[1])
colnames(Ggibbs) <- paste(species[site.all[,'specindex']],
                    SITE[site.all[,'SITE']],sep='.')
Ggibbs2 <- Ggibbs

Jgibbs  <- matrix(0,nt,dim(site.all)[1])
colnames(Jgibbs) <- paste(species[site.all[,'specindex']],
                    SITE[site.all[,'SITE']],sep='.')
Jgibbs2 <- Jgibbs

Wgibbs  <- matrix(0,nt,dim(site.all)[1])
colnames(Wgibbs) <- paste(species[site.all[,'specindex']],
                    SITE[site.all[,'SITE']],sep='.')
Wgibbs2 <- Wgibbs

ul      <- rep(0,nspec)
um      <- ul
ub      <- ul
ubs     <- ul
ug      <- ul
uk      <- ul  


