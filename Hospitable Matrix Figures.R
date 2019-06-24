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


###Persistence IN----

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



###Persistence Out----

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

