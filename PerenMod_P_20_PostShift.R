#100 years of population stabilization after climate shift
#Output:Map of stable allele distribution, text file of population numbers, saved workspace

load("PerenMod_P_20_Shift.RData")

#allow 100 yrs to re-stabilize
   TT3 <- 1:100  

step <- rep(0,length(KK))  #How far from ideal?
for(k in KK){
		step[k] <- abs(Ideal-kq[k])
}

for(t in TT3){ 
 #Plant area 
   PA <- matrix(0,length(KK),length(JJ))  
   for(k in 1:length(KK)){                #for each cell, for each size class...
   for(j in 1:length(JJ)){
   PA[k,j] <- sum(N[k,j,]*AJJ[j])    #number in size class x average basal area/class
   }}
	
   Comp.A <- matrix(0,length(KK),length(JJ))  # competitive area
   #for each size class...
   for(j in 1:(length(JJ)-1)){
   Comp.A[,j] <- apply(PA[,j:length(JJ)],1,sum)
   }
	
	for(k in KK){
	   #germination and survival 
			
			if(t > 1){                  #no seeds to germinate, initially
			germ <- max.germ*(TGclim[step[k]+1]-(0.003*sum(PA[k,])))
						
			#Seedling germination
			NewSdl <- 0
			if(germ > 0) NewSdl <- rbinom(1,sd.cell[k],germ)		
			N[k,1,] <- N[k,1,]+NewSdl
			} #end timestep loop
			
			#Survival
			for(j in JJ){
			surv <- 0	
				surv <- max.surv[j]*(SFclim[step[k]+1]-(CA.eff[j]*Comp.A[k,j]))
			if(surv<0) surv <- 0 
			N[k,j,]<- rbinom(1,N[k,j,],surv)
   			} #end size loop
	} #end cell loop
	
   #no seed bank - reset to 0
   sd.cell <- rep(0,length(KK))    
	
	#Producing seed 
	Ov.num <- rep(0,length(KK))  
	
	for(k in KK){
			fec <- rep(0,4)
	  if(sum(N[k,2:4,])>0){	
		#fecundity
		fec <- max.fec*SFclim[step[k]+1]
		Ov.num[k] <- sum(N[k,2:4,]*fec[2:4])
				} # end "if adults"
	}
     
  #seed dispersal from each cell to all others
   for (k in KK){
   	 if(sum(N[k,2:4,])>0){
   	    sd.cell <- sd.cell + Ov.num[k]*SD.mat[k,]
   	 }
  }  
  sd.cell <- round(sd.cell) #round numbers of seeds to whole numbers
    
   for(k in KK){
		for(i in II){               #transition to next size class for each genotype 
			for(j in 1:3){           #..and each size
			trans <- 0
			  trans <- max.trans[j]*(TGclim[step[k]+1]-(CA.eff[j]*Comp.A[k,j]))
			  if(trans<0) trans <- 0
			  			
			    Tr <- rbinom(1,N[k,j,i],trans); 
			    N[k,j,i] <- N[k,j,i]-Tr; N[k,j+1,i] <- N[k,j+1,i]+Tr
			} #end stage loop
		} #end genotype loop
	} #end cell loop

	print(t)
 }   #end setup cycle
 
#export as text file
Pop.i <- N

write.table(Pop.i,"PostShiftPop_Peren_P_20.txt")  
#write.table(Pop.i,"PostShiftPop_Peren_P_80.txt")    
 

   ##where are the genotypes (adults)?  
   pop.state <- rep(0,length(KK))
   occ.pshift <- rep(0,length(KK))
      
   for(k in KK){
   	if(sum(N[k,2:4,II])>0){
   		 occ.pshift[k] <- 1
   		 pop.state[k] <-  1 
   }}	
      		 
   mat.pop.state <- matrix(pop.state,ncol=Plot.width,nrow=Plot.length,byrow=T)
   jpeg(filename="PostShiftPop_Peren_P_alleles_20.jpg",width=9,height=7, units='in', res=500)
#jpeg(filename="PostShiftPop_Peren_P_alleles_80.jpg",width=9,height=7, units='in', res=500)
par(mfrow=c(2,2),mar=c(4,4,3,2))  
 image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=mat.pop.state,col=c("white","green4"),breaks=c(0,0.7,1),xlab="Distance",ylab="Distance",main="A1") 
 legend(50,950,legend=c("absent","present"), fill=c("white","green4") ,cex=0.7) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=clim.state ,col=c("yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3","violetred3","violetred4"),breaks=c(2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5),xlab="Distance",ylab="Distance",main="Climate") 
legend(50,950,legend=c("3","4","5","6","7","8","9","10","11","12"),fill=c("yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3","violetred3","violetred4"),cex=0.7)    
dev.off() 

#Change in occupancy
change.occ.ps <- sum(occ.pshift)/sum(occ)
percent.occ.ps <- sum(occ.pshift)/5000 

save.image(file="PerenMod_P_20_PostShift.RData")   
#save.image(file="PerenMod_P_80_PostShift.RData") 



