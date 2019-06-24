
### Set up Parameters----

#The variance for the randomized dispersal
sd <- .03
#The growth parameters (different for each species in the fish, same and lower in the water)
#These are just the default
Growth1 <- 1.4
Growth2 <- 1.45
Wa=1/2000 # set the carrying capacity in the water
a <- 1/200 # Beverton-holt alpha for all species 
init <- 60 # Initial colonization density
npatch<-seq(2,20,1) #list of patch numbers
exrate<-seq(.01,.1,.005) #list of patch disturbance rates
runs=100
time=2500
remain=0.01 #proportion of dispersers that remains in the water after dispersal event

### Set up backup functions -----

##Birth Function for all the patches at once
b <- function(R,current, patches, species, a) {
  birth <- matrix(NA, nrow=species, ncol=patches)
  for (p in 1:patches) {
    birth[,p] <- R[,p]*current[,p]*bevHoltC(current[,p],a)
  }
  return(birth)
}

## Competition within a given fish
bevHoltC <- function(N,a) {
  return( 1 / ( 1 + a*(sum(N))))
}


##Here I'm going to set up 2 new function -- one for direct one indirect
##These are going to be able to just set the traits in versus out

##TwoTraitDir Fish to fish or global   -----
TwoTraitDir <- function(patches,time,extinctionRate,outd,ind,samp=TRUE){
  
  #First Create the Data Array
  dat <- array(NA, c(2, patches, time))
  dat[,,1] <- init
  
  #Then Create The R matrix
  R <- matrix(NA, nrow=2, ncol=patches)
  R[1, ] <- Growth1
  R[2, ] <- Growth2
  
  
  #Now Set up the random extinction
  randommat<-matrix(NA,ncol=(patches),nrow=time)
  
  for(i in 1:(patches)){
    randommat[,i]<-rpois(time,extinctionRate)
  }
  
  
  #initialize, is there an extinction in the first time step 
  for (p in 1:(patches)) {
    if(randommat[1,p]==1) {
      dat[,p,1] <- rep(0,2)
    }
  }
  
  
  #Now start the cycles  
  for (t in 1:(time-1)){
    
    # birth
    
    birth <- b(R, dat[,,t],patches, 2, a)
    # # dispersal
    
    emigrate <- outd*birth[,c(1:patches)]
    stay <- birth[,c(1:patches)] - emigrate
    immigrate <- ind*emigrate
    #print(immigrate)
    
    # Updated Migration #---------------
    
    
    
    #This way does it so that the stuff that is direct fish to fish 
    
    if(samp==TRUE) dat[,,t+1] <- stay + immigrate[,sample(seq(1,patches))]
    else dat[,,t+1]<-stay + rowMeans(immigrate)
    #print(rowMeans(migrate))
    
    
    
    # extinction
    for (p in 1:patches) {
      if(randommat[t+1,p]==1) {
        dat[,p,t+1] <- rep(0,2)
      }
    }
  }
  
  return(dat) 
  
}
##TwoTraitInd Indirect through matrix ----
TwoTraitInd <- function(patches,time,extinctionRate,outd,ind,CC=FALSE){
  
  #First Create the Data Array
  dat <- array(NA, c(2, patches+1, time))
  dat[,,1] <- init
  
  #Then Create The R matrix
  R <- matrix(NA, nrow=2, ncol=patches+1)
  R[1, ] <- Growth1
  R[2, ] <- Growth2
  R[,patches+1] <- GrowthWater
  #print(R)
  #Now Set up the random extinction
  randommat<-matrix(NA,ncol=(patches+1),nrow=time)
  
  for(i in 1:(patches+1)){
    randommat[,i]<-rpois(time,extinctionRate)
  }
  
  
  #initialize, is there an extinction in the first time step 
  for (p in 1:(patches)) {
    if(randommat[1,p]==1) {
      dat[,p,1] <- rep(0,2)
    }
  }
  
  
  #Now start the cycles  
  for (t in 1:(time-1)){
    
    # birth
    
    birth <- b(R, dat[,,t],patches+1, 2, a)
    
    #print(birth)
    if(CC==TRUE) {birth[,patches+1] <- R[,patches+1]*dat[,patches+1,t]*bevHoltC(dat[,patches+1,t],Wa)}
    
    #print(birth)
    # # dispersal
    
    migrate <- outd*birth[,c(1:patches)]
    stay <- birth[,c(1:patches)] - migrate
    #print(migrate)
    
    # Updated Migration #---------------
    totalmig <- rowSums(migrate)
    waterpool <- birth[,patches+1] + rowSums(migrate)
    percentin <- sum(totalmig)/sum(waterpool)
    
    #print(percentin)
    ## In addition, I think we are going to need to remove a constant
    ## amount from the pool, which changes some of the way things happen
    
    
    dat[,c(1:patches),t+1] <- stay + ind*waterpool/patches
    print(ind*waterpool/patches)
    
    
    dat[,patches+1,t+1] <- waterpool-ind*waterpool
    
    # extinction
    for (p in 1:patches) {
      if(randommat[t+1,p]==1) {
        dat[,p,t+1] <- rep(0,2)
      }
    }
  }
  
  return(dat) 
  
}
##TwoTraitDirPool Fish to fish or global with a matrix-like patch  -----
TwoTraitDirPool <- function(patches,time,extinctionRate,outd,ind,wgrow,samp=TRUE){
  
  #First Create the Data Array
  dat <- array(NA, c(2, patches+1, time))
  dat[,,1] <- init
  
  #Then Create The R matrix
  R <- matrix(NA, nrow=2, ncol=(patches+1))
  R[1, ] <- Growth1
  R[2, ] <- Growth2
  R[,patches+1]<-wgrow
  # print(dim(R))
  #Now Set up the random extinction
  randommat<-matrix(NA,ncol=(patches),nrow=time)
  
  for(i in 1:(patches)){
    randommat[,i]<-rpois(time,extinctionRate)
  }
  
  
  #initialize, is there an extinction in the first time step 
  for (p in 1:(patches)) {
    if(randommat[1,p]==1) {
      dat[,p,1] <- rep(0,2)
    }
  }
  
  
  #Now start the cycles  
  for (t in 1:(time-1)){
    
    # birth
    
    birth <- b(R, dat[,,t],patches+1, 2, a)
    # # dispersal
    #print(birth)
    emigrate <- outd*birth
    #print(emigrate)
    stay <- birth - emigrate
    immigrate <- ind*emigrate
    
    #print(dim(immigrate))
    #print(migrate)
    
    # Updated Migration #---------------
    
    
    
    #This way does it so that the stuff that is direct fish to fish 
    
    if(samp==TRUE) dat[,,t+1] <- stay + immigrate[,sample(seq(1,patches+1))]
    else dat[,,t+1]<-stay + rowMeans(immigrate)
    #print(rowMeans(migrate))
    
    
    
    # extinction
    for (p in 1:patches) {
      if(randommat[t+1,p]==1) {
        dat[,p,t+1] <- rep(0,2)
      }
    }
  }
  
  return(dat) 
  
}


