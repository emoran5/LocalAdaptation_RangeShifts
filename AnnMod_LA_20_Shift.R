#Climate will shift at 20 m/yr for 20 years.  Climates 1 & 2 (marginally suitable) will entirely disappear.  Lowest 2 bands will be filled with climate 11 (very marginally suitable) or 12 (not suitable for species) 
#Output:Map of stable allele distribution, text file of population numbers, saved workspace

load("AnnMod_LA_20_Stabil.RData")

#Duration of climate shift
   TT2 <- 1:20

for(t in TT2){   
#climate shift
	if(t <= 10) kq <- c(rep(11,Plot.width),kq[1:(CellNum-Plot.width)])
	if(t > 10) kq <- c(rep(12,Plot.width),kq[1:(CellNum-Plot.width)])
		
	step <- matrix(0,length(KK),6)  #How far from ideal?
for(k in KK){
	for(i in II){
		step[k,i] <- abs(Ideal[i]-kq[k])
	}
}

N[,2,] <- 0  #no adults survive to next year
	   	
	for(k in KK){
		for(i in II){               #germination and survival for each genotype (add competition)
		
	##Germination - affected by climate	
				germ <- 0	
				germ<- max.germ*TGclim[step[k,i]+1]
            
			NewSdl <- 0
			if(germ > 0) NewSdl <- rbinom(1,N[k,1,i],germ)	
			
			N[k,1,i] <- N[k,1,i]-NewSdl
	
	##Climate impact on survival
	    surv <- 0
	    
	if (i ==1) surv<- SFclim.1.6[step[k,i]+1]
	if (i ==6) surv<- SFclim.1.6[step[k,i]+1]
	if (i>1 & i<6) surv <- surv<- SFclim.2.5[step[k,i]+1]
	
	NewSdl2 <- rbinom(1,NewSdl,surv)
			
			##Survival of seeds to next year
			N[k,1,i] <- 0 #if no seedbank
			
			##Adult placeholder
			N[k,2,i] <- N[k,2,i]+NewSdl2
			
		} #end genotype loop
	
	##Overall density effects on survival	
		surv <- 0
		Dens.adj <- 0
		
		Sdl.dens <- sum(N[k,2,])
		if(7000/Sdl.dens >= 1) surv <- 1
		if(7000/Sdl.dens < 1) surv <- 7000/Sdl.dens

		if(surv <1) for (i in II){
			N[k,2,i]<- rbinom(1,N[k,2,i],surv)
		}
		
	} #end cell loop    
	
	#Producing seed and pollen
	Ov.num <- matrix(0,length(KK),3);   #Number of A1, A2, or A3 ovules in cell K
	Pol.num <- matrix(0,length(KK),3);  #Number of A1, A2, or A3 pollen in cell K
	
	for(k in KK){
	  if(sum(N[k,2,])>0){	
		#fecundityfor each genotype
		fec <- rep(0,6); pol <- rep(0,6)
		
		for(i in II){   
			if(i==1){
				fec[i] <- max.fec*SFclim.1.6[step[k,i]+1]
				pol[i] <- max.pol*SFclim.1.6[step[k,i]+1]		
			}            
			if(i==6){
				fec[i] <- max.fec*SFclim.1.6[step[k,i]+1]
				pol[i] <- max.pol*SFclim.1.6[step[k,i]+1]		
			}  
			if(i>1 & i<6){
				fec[i] <- max.fec*SFclim.2.5[step[k,i]+1]
				pol[i] <- max.pol*SFclim.2.5[step[k,i]+1]		
			}  
		} # end "if adults"
		   
			Ov.num[k,1] <- (sum(N[k,2,1]*fec[1])+0.5*sum(N[k,2,2]*fec[2])+0.5*sum(N[k,2,4]*fec[4]))  #A1
			Ov.num[k,2] <- (sum(N[k,2,3]*fec[3])+0.5*sum(N[k,2,2]*fec[2])+0.5*sum(N[k,2,5]*fec[5]))  #A2
			Ov.num[k,3] <- (sum(N[k,2,6]*fec[6])+0.5*sum(N[k,2,4]*fec[4])+0.5*sum(N[k,2,5]*fec[5]))  #A3
			
			Pol.num[k,1] <- (sum(N[k,2,1]*pol[1])+0.5*sum(N[k,2,2]*pol[2])+0.5*sum(N[k,2,4]*pol[4]))  #A1
			Pol.num[k,2] <- (sum(N[k,2,3]*pol[3])+0.5*sum(N[k,2,2]*pol[2])+0.5*sum(N[k,2,5]*pol[5]))  #A2
			Pol.num[k,3] <- (sum(N[k,2,6]*pol[6])+0.5*sum(N[k,2,4]*pol[4])+0.5*sum(N[k,2,5]*pol[5]))  #A3
			
    } #end "if adults present" loop
	} #end k loop
   
   #Pollen dispersal
   #columns - pollen alleles
   Pol.rec <- matrix(0,length(KK),3)  
   for (k in KK){
   	  if(sum(N[k,2,])>=1){ 
   	    Pol.rec[,1] <- Pol.rec[,1] + Pol.num[k,1]*Pol.mat[k,] #number of a1 pollen grains dispersed to each  other square
   	    Pol.rec[,2] <- Pol.rec[,2] + Pol.num[k,2]*Pol.mat[k,] #number of a1 pollen grains dispersed to each  other square
   	    Pol.rec[,3] <- Pol.rec[,3] + Pol.num[k,3]*Pol.mat[k,] #number of a1 pollen grains dispersed to each  other square
   	  }
   }
   
   #proportion of received pollen of each genotype
   Pol.prop.a1 <- rep(0,length(KK))
   Pol.prop.a2 <- rep(0,length(KK))
   Pol.prop.a3 <- rep(0,length(KK))
   
   for(i in KK){
   Pol.prop.a1[i] <- Pol.rec[i,1]/sum(Pol.rec[i,1:3])
   Pol.prop.a2[i] <- Pol.rec[i,2]/sum(Pol.rec[i,1:3])
   Pol.prop.a3[i] <- Pol.rec[i,3]/sum(Pol.rec[i,1:3])   	
   }

   #seed genotype numbers
   S.fert <- matrix(0,length(KK),6)
     S.fert[,1]<- Ov.num[,1]*Pol.prop.a1 
     S.fert[,2]<- Ov.num[,2]*Pol.prop.a1 + Ov.num[,1]*Pol.prop.a2
     S.fert[,3]<- Ov.num[,2]*Pol.prop.a2 
     S.fert[,4]<- Ov.num[,3]*Pol.prop.a1 + Ov.num[,1]*Pol.prop.a3
     S.fert[,5]<- Ov.num[,3]*Pol.prop.a2 + Ov.num[,2]*Pol.prop.a3
     S.fert[,6]<- Ov.num[,3]*Pol.prop.a3 
     
  #seed dispersal from each cell to all others
   for (k in KK){
   	 if(sum(N[k,2,])>0){
   	    N[,1,1] <- N[,1,1] + S.fert[k,1]*SD.mat[k,]
   	    N[,1,2] <- N[,1,2] + S.fert[k,2]*SD.mat[k,]
   	    N[,1,3] <- N[,1,3] + S.fert[k,3]*SD.mat[k,]
   	    N[,1,4] <- N[,1,4] + S.fert[k,4]*SD.mat[k,]
   	    N[,1,5] <- N[,1,5] + S.fert[k,5]*SD.mat[k,]
   	    N[,1,6] <- N[,1,6] + S.fert[k,6]*SD.mat[k,]
   	 }
  }
  
  N[,1,] <- round(N[,1,]) #round numbers of seeds to whole numbers

	print(t)
		
}	#end cycle	

#export as text file
Pop.i <- cbind(N[,,1],N[,,2],N[,,3],N[,,4],N[,,5],N[,,6])
#write.table(Pop.i,"ShiftPop_Annual_LA_1.txt")  
write.table(Pop.i,"ShiftPop_Annual_LA_20.txt")   
#write.table(Pop.i,"ShiftPop_Annual_LA_50.txt") 

 ##where are the genotypes (adults)?  
   pop.state <- rep(0,length(KK))
   prop1 <- rep(0,length(KK)) #proportion of A1
   prop2 <- rep(0,length(KK)) #proportion of A2
   prop3 <- rep(0,length(KK)) #proportion of A3
   occ.shift <- rep(0,length(KK))		
   		
   for(k in KK){
   	if(sum(N[k,2,II])>0){
   		 occ.shift[k] <- 1
   		 prop1[k] <-  (sum(0.5*N[k,2,2])+sum(0.5*N[k,2,4])+sum(N[k,2,1]))/sum(N[k,2,II]) 
   		 prop2[k] <-  (sum(0.5*N[k,2,5])+sum(0.5*N[k,2,2])+sum(N[k,2,3]))/sum(N[k,2,II]) 
   		 prop3[k] <-  (sum(0.5*N[k,2,4])+sum(0.5*N[k,2,5])+sum(N[k,2,6]))/sum(N[k,2,II]) 
   		  	
   	if(prop1[k] >= 0.75 & prop2[k] <= 0.25) pop.state[k] <- 1  #>=75% A1 alleles, <= 25% A2 alleles
   	if(prop1[k] >= 0.55 & prop1[k] <= 0.75 & prop2[k] <= 0.3 & prop3[k] <= 0.15) pop.state[k] <- 2  #>=55% A1 alleles, <= 30% A2 alleles, <= 15% A3 alleles
   	if(prop1[k] >= 0.35 & prop1[k] <= 0.55 & prop2[k] >= 0.35 & prop3[k] <= 0.30) pop.state[k] <- 3  #>=35% A1 alleles, >= 35% A2 alleles, <= 30% A3 alleles
    if(prop1[k] <= 0.3 & prop2[k] >= 0.4 & prop3[k] <= 0.30) pop.state[k] <- 4  #<=30% A1 alleles, >= 40% A2 alleles, <= 30% A3 alleles 	
   	if(prop1[k] <= 0.3 & prop2[k] >= 0.35 & prop3[k] <= 0.55 & prop3[k] >= 0.35) pop.state[k] <- 5
   	  #<= 30% A1 alleles, >= 35% A2 alleles, >=35%  A3 alleles
    if(prop1[k] <= 0.15 & prop2[k] <= 0.3 & prop3[k] <= 0.75 & prop3[k] >= 0.55) pop.state[k] <- 6  #>=55% A3 alleles, <= 30% A2 alleles, <= 15% A1 alleles
   	if(prop2[k] <= 0.25 & prop3[k] >= 0.75) pop.state[k] <- 7  #>= 70% A3 alleles, <= 30% A2 alleles
   	}}
   		 
   mat.pop.state <- matrix(pop.state,ncol=Plot.width,nrow=Plot.length,byrow=T)
   prop1.mat  <- matrix(prop1,ncol=Plot.width,nrow=Plot.length,byrow=T)
   prop2.mat  <- matrix(prop2,ncol=Plot.width,nrow=Plot.length,byrow=T)
   prop3.mat  <- matrix(prop3,ncol=Plot.width,nrow=Plot.length,byrow=T) 
   
 clim.state <- matrix(kq,ncol=Plot.width,nrow=Plot.length,byrow=T)   

#jpeg(filename="ShiftPop_Annual_LA_alleles_1.jpg",width=9,height=7, units='in', res=500)
jpeg(filename="ShiftPop_Annual_LA_alleles_20.jpg",width=9,height=7, units='in', res=500)
#jpeg(filename="ShiftPop_Annual_LA_alleles_50.jpg",width=9,height=7, units='in', res=500)
par(mfrow=c(2,2),mar=c(4,4,3,2))  
 image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop1.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A1") 
 legend(50,950,legend=c("<0.1","0.1-0.3","0.3-0.5","0.5-0.7",">0.7"), fill=c("white","yellow","gold2","green2","green4") ,cex=0.7) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop2.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A2") 
 
image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop3.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A3") 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=clim.state ,col=c("yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3","violetred3","violetred4"),breaks=c(2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5),xlab="Distance",ylab="Distance",main="Climate") 
legend(50,950,legend=c("3","4","5","6","7","8","9","10","11","12"), fill=c("yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3","violetred3","violetred4"),cex=0.7)  
dev.off()  

#Change in occupancy
change.occ <- sum(occ.shift)/sum(occ)
percent.occ <- sum(occ.shift)/5000 		

#save.image(file="AnnMod_LA_1_Shift.RData")   
save.image(file="AnnMod_LA_20_Shift.RData") 
#save.image(file="AnnMod_LA_50_Shift.RData")  