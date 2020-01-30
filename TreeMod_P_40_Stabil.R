#200 years of population stabilization
#Output:Map of stable allele distribution, text file of population numbers, saved workspace

load("TreeMod_P_Setup.RData")

 #Setup initial conditions
   TT1 <- 1:200
   
   #Which seed dispersal?
SD.mat <- SD.mat.40
#SD.mat <- SD.mat.80


N <- N.i  #initial population

   #Dispersed seeds in the beginning of the simulation - none
   sd.cell <- rep(0,length(KK))     

step <- rep(0,length(KK))  #How far from ideal?
for(k in KK){
		step[k] <- abs(Ideal-kq[k])
}

for(t in TT1){  
	  #Basal area 
   BA <- matrix(0,length(KK),length(JJ))  
   for(k in 1:length(KK)){                #for each cell, for each size class...
   for(j in 1:length(JJ)){
   BA[k,j] <- sum(N[k,j,]*BAJJ[j])    #number in size class x average basal area/class
   }}
	
   Comp.BA <- matrix(0,length(KK),length(JJ))  # competitive basal area
   #for each size class...
   for(j in 1:(length(JJ)-1)){
   Comp.BA[,j] <- apply(BA[,j:length(JJ)],1,sum)
   }
	
	for(k in KK){
		#germination and survival 
			
			if(t > 1){              #no seeds to germinate, initially
			germ <- max.germ*(TGclim[step[k]+1]-(0.08*sum(BA[k,])))
						
			#Seedling germination
			NewSdl <- 0
			if(germ > 0) NewSdl <- rbinom(1,sd.cell[k],germ)		
			N[k,1,] <- N[k,1,]+NewSdl
			} #end timestep loop
			
			#Survival
			for(j in JJ){
			surv <- 0	
				surv <- max.surv[j]*(SFclim[step[k]+1]-(BA.eff[j]*Comp.BA[k,j]))
				if(surv<0) surv <- 0 
			N[k,j,]<- rbinom(1,N[k,j,],surv)
   			} #end size loop
	} #end cell loop
	
   #no seed bank - reset to 0
   sd.cell <- rep(0,length(KK))     
	
	#Producing seed 
	Ov.num <- rep(0,length(KK))   
		
		for(k in KK){
			fec <- rep(0,5)
	  if(sum(N[k,3:5,])>0){	
		#fecundity
		fec <- max.fec*SFclim[step[k]+1]
		Ov.num[k] <- sum(N[k,3:5,]*fec[3:5])
				} # end "if adults"
	}
   
  #seed dispersal from each cell to all others
   for (k in KK){
   	 if(sum(N[k,3:5,])>0){
   	    sd.cell <- sd.cell + Ov.num[k]*SD.mat[k,]
   	 }
  }  
  sd.cell <- round(sd.cell) #round numbers of seeds to whole numbers
  
   for(k in KK){
		#transition to next size class 
			for(j in 1:4){           #..and each size
			trans <- 0
			  trans <- max.trans[j]*(TGclim[step[k]+1]-(BA.eff[j]*Comp.BA[k,j]))
			  if(trans<0) trans <- 0
			  			
			    Tr <- rbinom(1,N[k,j,],trans); 
			    N[k,j,] <- N[k,j,]-Tr; N[k,j+1,] <- N[k,j+1,]+Tr
			} #end stage loop
	} #end cell loop

	print(t)
 }   #end setup cycle
 
#export as text file
Pop.i <- N
write.table(Pop.i,"StablePop_Tree_P_40.txt")  
#write.table(Pop.i,"StablePop_Tree_P_80.txt")    

   ##where are the genotypes (adults)?  
   pop.state <- rep(0,length(KK))
   occ <- rep(0,length(KK))	
      
   for(k in KK){
   	if(sum(N[k,3:5,II])>0){
   		 occ[k] <- 1
   		 pop.state[k] <-  1 
   }}		 
      		 
mat.pop.state <- matrix(pop.state,ncol=Plot.width,nrow=Plot.length,byrow=T) 
   
   jpeg(filename="StablePop_Tree_P_40.jpg",width=9,height=7, units='in', res=500)
#jpeg(filename="StablePop_Tree_P_80.jpg",width=9,height=7, units='in', res=500)
par(mfrow=c(2,2),mar=c(4,4,3,2))  
image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=mat.pop.state,col=c("white","green4"),breaks=c(0,0.7,1),xlab="Distance",ylab="Distance",main="A1") 
 legend(50,950,legend=c("absent","present"), fill=c("white","green4") ,cex=0.7) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=clim.state ,col=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5),xlab="Distance",ylab="Distance",main="Climate") 
legend(50,950,legend=c("1","2","3","4","5","6","7","8","9","10"), fill=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),cex=0.7)  
dev.off()   

save.image(file="TreeMod_P_40_Stabil.RData")   
#save.image(file="TreeMod_P_80_Stabil.RData") 

