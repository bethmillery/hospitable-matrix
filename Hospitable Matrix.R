
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
    #print(migrate)
    
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
    #print(d*waterpool/patches)
    
    
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


# Figures in the paper -------
###
### Both Traits at once
### Figure 2
###

fold=6^.5 ## each trait gets Sqrt(6) fold difference since the effect is multiplicative
meanout=.01 #the mean rate of getting out of the fish Supplemental figure 4 alters this and meanin
meanin=.1 #The mean rate of getting into the fish from the matrix
varout=c(fold*meanout*2/(fold+1),meanout*2/(fold+1)) #Calculate the two rates if the trait difference is on the way out
varin=c(fold*meanin*2/(fold+1),meanin*2/(fold+1)) #Calculate the two rates if the trait difference is on the way in


###Direct----
# Create an arry to hold the data
denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))
## loop through patches, patch disturbance, and reps
for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      #use the direct dispersal function with Samp=True to create patch-to-patch dispersal
      #varout and varin set the dispersal specialist advantage to be in both traits
      outFF<-TwoTraitDir(npatch[np],time,exrate[ex],varout,varin,TRUE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}

##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")


###Global----

##create an array to hold the data
denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))
## loop through patches, patch disturbance, and reps
for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      
      #use the direct dispersal function with samp=False to create global dispersal
      outFF<-TwoTraitDir(npatch[np],time,exrate[ex],varout,varin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}

##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")



###Persistence----
#Set the Matrix Growth rate to 1 (no growth)
GrowthWater=c(1,1)
# create an array
denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))
## loop through patches, patch disturbance, and reps
for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      #Use the indirect dispersal function
      outFF<-TwoTraitInd(npatch[np],time,exrate[ex],varout,varin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}



##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")




###Growth ----
##Set growth in the water to be lower than in the patch but equal between species
GrowthWater=c(1.2,1.2)

denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      
      outFF<-TwoTraitInd(npatch[np],time,exrate[ex],varout,varin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}

##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")



#----
####
#### Fold = 6
#### Figure 4
fold=6 # now we will look at the full advantage in only one trait at a time (supplement S2-3 have set to 10 and 3)
meanout=.01 #the mean rate of getting out of the fish
meanin=.1 #The mean rate of getting into the fish from the matrix
varout=c(fold*meanout*2/(fold+1),meanout*2/(fold+1)) #Calculate the two rates if the trait difference is on the way out
varin=c(fold*meanin*2/(fold+1),meanin*2/(fold+1)) #Calculate the two rates if the trait difference is on the way in

###Direct OUT----


denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      
      #varout=variation in dout, meanin=no variation in din
      outFF<-TwoTraitDir(npatch[np],time,exrate[ex],varout,meanin,TRUE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
     
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")


###Direct IN----


denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      
      #meanout=no variation in dout, varin=variation in din
      outFF<-TwoTraitDir(npatch[np],time,exrate[ex],meanout,varin,TRUE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")



###Global OUT----


denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      
      #varout=variation in dout, meanin=no variation in din
      outFF<-TwoTraitDir(npatch[np],time,exrate[ex],varout,meanin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")


###Global IN----


denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      
      #meanout=no variation in dout, varin=variation in din
      outFF<-TwoTraitDir(npatch[np],time,exrate[ex],meanout,varin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")



###Low Growth OUT----

GrowthWater=c(1.2,1.2)

denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      #varout=variation in dout, meanin=no variation in din
      outFF<-TwoTraitInd(npatch[np],time,exrate[ex],varout,meanin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")


###Low Growth IN----

GrowthWater=c(1.2,1.2)

denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      #meanout=no variation in dout, varin=variation in din
      outFF<-TwoTraitInd(npatch[np],time,exrate[ex],meanout,varin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")


###Memory IN----

GrowthWater=c(1,1)

denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      #meanout=no variation in dout, varin=variation in din
      outFF<-TwoTraitInd(npatch[np],time,exrate[ex],meanout,varin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")



###Memory Out----

GrowthWater=c(1,1)

denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      #varout=variation in dout, meanin=no variation in din
      outFF<-TwoTraitInd(npatch[np],time,exrate[ex],varout,meanin,FALSE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")



###----
### Selection in the Patch vs. in the Matrix
### Figure 5----

fold=6^(1/2)
meanout=.01

meanin=.1

#fast growth in the matrix (adjust Wsp1 for slower matrix growth)

WSp2<-seq(1.4,1.6,.01)
WSp1<-seq(1.4,1.6,.01)

denseMat<-array(NA,dim=c(2,runs,length(WSp1),length(WSp2)))

for(f in 1:length(WSp1)){
  for(w2 in 1:length(WSp2)){
    for(r in 1:runs){
      
      Growth2 <- WSp2[w2]
      GrowthWater=c(WSp1[f],1.4)
      
      #There is no variation in the din or dout traits
      #only 4 patches, adjust for additional patches
      outFF<-TwoTraitInd(4,time,.05,meanout,meanin,FALSE) 
      
      denseMat[1,r,f,w2]<-sum(outFF[1,1:4,c((time-50):time)])/(500)
      denseMat[2,r,f,w2]<-sum(outFF[2,1:4,c((time-50):time)])/(500)
      
      if(r==1){
        print(denseMat[,r,f,w2])
      }
      
    }
    print(c( WSp2[w2],WSp1[f]))
  }
}



##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=WSp2,y=WSp1,z=t(percentextinct[,1,]),main="Dispersal Specialist")
p1=filled.contour(x=WSp2,y=WSp1,z=t(percentextinct[,2,]),main="Patch Specialist")


###----
### Matrix-like patch
### Figure 6
###Direct with growth in the matrix-like patch----
fold=6^.5
meanout=.01
varout=c(fold*meanout*2/(fold+1),meanout*2/(fold+1))
meanin=.1
varin=c(fold*meanin*2/(fold+1),meanin*2/(fold+1))

denseMat<-array(NA,dim=c(2,runs,length(npatch),length(exrate)))

for(np in 1:length(npatch)){
  for(ex in 1:length(exrate)){
    for(r in 1:runs){
      #Growth in the matrix-like patch = c(1.2,1.2) (change to c(1,1) for persistence)
      #Use the function TwoTraitDirPool to implement
      outFF<-TwoTraitDirPool(npatch[np],time,exrate[ex],varout,varin,c(1.2,1.2),TRUE) 
      
      denseMat[1,r,np,ex]<-sum(outFF[1,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      denseMat[2,r,np,ex]<-sum(outFF[2,1:npatch[np],c((time-50):time)])/(50*npatch[np])
      
      if(r==1){
        
        print(denseMat[,r,np,ex])
      }
      
    }
    print(exrate[ex])
  }
  print(npatch[np])
  
}


##calculate the percent of times a species is found in the system
percentextinct=array(dim=c(length(npatch),2,length(exrate)))
#loop through the denseMat array
for(i in 1:length(npatch)){
  
  for(j in 1:length(exrate)){
    
    percentextinct[i,1,j]<-sum(denseMat[1,,i,j]>1,na.rm=TRUE)
    percentextinct[i,2,j]<-sum(denseMat[2,,i,j]>1,na.rm=TRUE)
    
  }
}

##quick plot
par(mar=c(3,3,3,4))
p1=filled.contour(x=npatch,y=exrate,z=percentextinct[,1,],main="Dispersal Specialist")
p2=filled.contour(x=npatch,y=exrate,z=percentextinct[,2,],main="Patch Specialist")

