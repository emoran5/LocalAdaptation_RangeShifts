#20 years of population stabilization
#Output:Map of stable distribution, text file of population numbers, saved workspace

load("AnnMod_P_Setup.RData")

 #Setup initial conditions
   TT1 <- 1:20

#Which seed dispersal?
SD.mat <- SD.mat.20
#SD.mat <- SD.mat.1
#SD.mat <- SD.mat.50   

N <- N.i  #initial population   

step <- rep(0,length(KK))  #How far from ideal?
for(k in KK){
		step[k] <- abs(Ideal-kq[k])
}

for(t in TT1){  
	 
	   	N[,2,] <- 0  #no adults survive to next year
	   	
	for(k in KK){
#germination and survival 	##Germination - affected by climate	
				germ <- 0	
				germ<- max.germ*TGclim[step[k]+1]
            
			NewSdl <- 0
			if(germ > 0) NewSdl <- rbinom(1,N[k,1,],germ)	
			
			N[k,1,] <- N[k,1,]-NewSdl
	
	##Climate impact on survival
	    surv <- 0
	    
	surv<- SFclim[step[k]+1]
	
	NewSdl2 <- rbinom(1,NewSdl,surv)
			
			##Survival of seeds to next year
			N[k,1,] <- 0 #if no seedbank
			
			##Adult placeholder
			N[k,2,] <- N[k,2,]+NewSdl2
			
	
	##Overall density effects on survival	
		surv <- 0
		Dens.adj <- 0
		
		Sdl.dens <- sum(N[k,2,])
		if(7000/Sdl.dens >= 1) surv <- 1
		if(7000/Sdl.dens < 1) surv <- 7000/Sdl.dens

		if(surv <1) N[k,2,i]<- rbinom(1,N[k,2,i],surv)
		
	} #end cell loop    
	
	#Producing seed 
	Ov.num <- rep(0,length(KK));
	
	for(k in KK){
	  if(sum(N[k,2,])>0){	
		#fecundity
		Ov.num[k] <- max.fec*SFclim[step[k]+1]
				} # end "if adults"
	}
     
  #seed dispersal from each cell to all others
   for (k in KK){
   	 if(sum(N[k,2,])>0){
   	    N[,1,] <- N[,1,] + Ov.num[k]*SD.mat[k,]
   	 }
  }
  
  N[,1,] <- round(N[,1,]) #round numbers of seeds to whole numbers

	print(t)
 }   #end setup cycle
 
#export as text file
Pop.i <- N
#write.table(Pop.i,"StablePop_Annual_P_1.txt")  
write.table(Pop.i,"StablePop_Annual_P_20.txt")   
#write.table(Pop.i,"StablePop_Annual_P_50.txt") 
 
   ##where are the genotypes (adults))?  
   pop.state <- rep(0,length(KK))
   occ <- rep(0,length(KK))		
   		
   for(k in KK){
   	if(sum(N[k,2,II])>0){
   		 occ[k] <- 1
   		 pop.state[k] <-  1 
   	}}
   		 
mat.pop.state <- matrix(pop.state,ncol=Plot.width,nrow=Plot.length,byrow=T)   
 
#jpeg(filename="StablePop_Annual_P_1.jpg",width=9,height=7, units='in', res=500)
jpeg(filename="StablePop_Annual_P_20.jpg",width=9,height=7, units='in', res=500)
#jpeg(filename="StablePop_Annual_P_50.jpg",width=9,height=7, units='in', res=500)
par(mfrow=c(2,2),mar=c(4,4,3,2))  
image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=mat.pop.state,col=c("white","green4"),breaks=c(0,0.7,1),xlab="Distance",ylab="Distance",main="A1") 
 legend(50,950,legend=c("absent","present"), fill=c("white","green4") ,cex=0.7) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=clim.state ,col=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5),xlab="Distance",ylab="Distance",main="Climate") 
legend(50,950,legend=c("1","2","3","4","5","6","7","8","9","10"), fill=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),cex=0.7) 
dev.off()

#save.image(file="AnnMod_P_1_Stabil.RData")   
save.image(file="AnnMod_P_20_Stabil.RData") 
#save.image(file="AnnMod_P_50_Stabil.RData")  