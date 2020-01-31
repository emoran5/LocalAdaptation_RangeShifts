#50 x 100 20mx20m plot landscape w/10 climates.  
# Columns 1-10 climate 10, col 11-20 climate 9, col 21-30 climate 8; col 31-40 climate 7; col 41-50 climate 6; col 51-60 climate 5; col 61-70 climate 4; col 71-80 climate 3; col 81-90 climate 2; and col 91-100 climate 1.

#For genotype-specific simulation, varying proportions of each genotype by climate band.

#Output: Map of initial allele distribution, text file of population numbers, dispersal kernel graphs, dispersal matrices saved in R file 

##############Setup plot
Plot.width <- 50
Plot.length <- 100
C.width <- Plot.length/10
CellNum <- (Plot.width*Plot.length)

 #Habitat cells
   KK <- 1:CellNum
   #habitat cell quality - climate, given by 1 (cold) to 10 (hot)
kq <- c(rep(10,Plot.width*C.width),rep(9,Plot.width*C.width),rep(8,Plot.width*C.width),rep(7,Plot.width*C.width),rep(6,Plot.width*C.width),rep(5,Plot.width*C.width),rep(4,Plot.width*C.width),rep(3,Plot.width*C.width),rep(2,Plot.width*C.width),rep(1,Plot.width*C.width)) 
   
   #Cell indexing - form into rectangle
   index <- cbind(row=rep(1:Plot.length,each=Plot.width),column=rep(1:Plot.width,Plot.length)) #xy coord
   
 clim.state <- matrix(kq,ncol=Plot.width,nrow=Plot.length,byrow=T)
    

########Annual population
 #Two size classes - seeds and adults
   JJ <- 1:2

   #Genotypes 
   II <- 1:6   
   #A1A1; A1A2; A2A2; A3A1; A2A3; A3A3
   
      #Initial population structure
   n.entry <- length(KK)*length(II)*length(JJ) 
   #habitat cells (5000) x genotypes (6) x size classes (5)
   N.i <- rep(0,n.entry)
   dim(N.i) <- c(length(KK),length(JJ),length(II)) #5000 x 5 matrix, set of 6
   
   ##6 genotypes, and where they start out
   for(k in KK){
   #Climate 1...10
   	if(kq[k]==10){
   N.i[k,1,] <- c(0,0,0,0,100,300)   #seeds - 400
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
   	}
   	if(kq[k]==9){
   N.i[k,1,] <- c(0,0,20,20,260,500)   #seeds - 800
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
   	}
   	if(kq[k]==8){
   N.i[k,1,] <- c(0,100,150,150,600,1000)   #seeds - 2000
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
      	} 
   	if(kq[k]==7){
   N.i[k,1,] <- c(0,200,250,250,1000,300)    #seeds - 2000
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
      	}
   	if(kq[k]==6){
   N.i[k,1,] <- c(100,400,450,450,500,100)    #seeds - 2000
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
   	}
   	if(kq[k]==5){
   N.i[k,1,] <- c(300,1000,250,250,200,0)    #seeds - 2000
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
   	}
   	if(kq[k]==4){
   N.i[k,1,] <- c(1000,600,150,150,100,0)    #seeds - 2000
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
      	}
   	if(kq[k]==3){
   N.i[k,1,] <- c(500,260,20,20,0,0)    #seeds - 800
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
   	}
   	if(kq[k]==2){
   N.i[k,1,] <- c(300,100,0,0,0,0)   #seeds - 400
   N.i[k,2,] <- c(0,0,0,0,0,0)  #adults
   	}
   }
   ##where are the genotypes (seeds)?  
   pop.state <- rep(0,length(KK))
   prop1 <- rep(0,length(KK)) #proportion of A1
   prop2 <- rep(0,length(KK)) #proportion of A2
   prop3 <- rep(0,length(KK)) #proportion of A3
   
   for(k in KK){
   	if(sum(N.i[k,1,II])>0){
   		 prop1[k] <-  (sum(0.5*N.i[k,1,2])+sum(0.5*N.i[k,1,4])+sum(N.i[k,1,1]))/sum(N.i[k,1,II]) 
   		 prop2[k] <-  (sum(0.5*N.i[k,1,5])+sum(0.5*N.i[k,1,2])+sum(N.i[k,1,3]))/sum(N.i[k,1,II]) 
   		 prop3[k] <-  (sum(0.5*N.i[k,1,4])+sum(0.5*N.i[k,1,5])+sum(N.i[k,1,6]))/sum(N.i[k,1,II]) 
   		  	
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
   
jpeg(filename="Initial_Annual_LA_alleles.jpg",width=9,height=7, units='in', res=500)
par(mfrow=c(2,2),mar=c(4,4,3,2))  
 image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop1.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A1") 
 legend(50,950,legend=c("<0.1","0.1-0.3","0.3-0.5","0.5-0.7",">0.7"), fill=c("white","yellow","gold2","green2","green4") ,cex=0.7) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop2.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A2") 
 
image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop3.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A3") 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=clim.state ,col=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5),xlab="Distance",ylab="Distance",main="Climate") 
legend(50,950,legend=c("1","2","3","4","5","6","7","8","9","10"), fill=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),cex=0.7)  
dev.off()    

#export as text file
Pop.i <- cbind(N.i[,,1],N.i[,,2],N.i[,,3],N.i[,,4],N.i[,,5],N.i[,,6])
write.table(Pop.i,"InitialPop_Annual_LA.txt")      


####### Vegetative processes
Ideal <- c(4,5,6,6,7,8)   #Ideal climate for each genotype 

max.fec <- 300  #maximum seed production for each size class #was 500
max.pol <- 8000  #maximum pollen production for each size class #was 10000
max.germ <- 0.3		#maximum germination 
#max.surv <- 0.05 #maximum survival for seeds in seedbank

SFclim.1.6 <- c(0.9,0.5,0.3,0.15,0,0,0,0,0) #survival/fecundity effect of climate, genotypes 1 & 6; ideal, 1-6 steps
SFclim.2.5 <- c(1,0.7,0.4,0.25,0.1,0,0,0,0) #survival/fecundity effect of climate, genotypes 2-5; ideal, 1-6 steps
#SFclim.1.6 <- c(0.9,0.75,0.55,0.4,0.3,0.2,0.1,0,0) #survival/fecundity effect of climate, genotypes 1 & 6; ideal, 1-6 steps
#SFclim.2.5 <- c(1,0.8,0.6,0.45,0.35,0.2,0.1,0,0) #survival/fecundity effect of climate, genotypes 2-5; ideal, 1-6 steps

TGclim <- c(1,0.7,0.5,0.2,0.1,0,0,0,0) #transition/germination effect of climate; ideal, 1 step-6 steps
#TGclim <- c(1,0.9,0.8,0.6,0.4,0.2,0.1,0,0) #transition/germination effect of climate; ideal, 1 step-6 steps

##########Dispersal
###2D-t seed dispersal
#mean dispersal = 1 m
u.1 <- (2/pi)^2 #0.41
#mean dispersal = 20 m
u.20 <- (20*2/pi)^2 #364.8
#mean dispersal = 50 m
u.50 <- (50*2/pi)^2 #1013.2

#pollen mean dispersal = 60 m
u.p <- (60*2/pi)^2 #1459

x1 <- seq(0,2010,by=0.1)
pd.1 <- 1/(pi*u.1*(1+((x1^2)/u.1))^2)
pd.20 <- 1/(pi*u.20*(1+((x1^2)/u.20))^2)
pd.50 <- 1/(pi*u.50*(1+((x1^2)/u.50))^2)
pd.p <- 1/(pi*u.p*(1+((x1^2)/u.p))^2)


jpeg(filename="DispersalKernels_annuals",width=6,height=6, units='in', res=500)
par(mfrow=c(2,1),mar=c(4,4,2,2))
plot(x1,pd.1,xlim=c(0,20),type='l',lwd=2,ylab="probability density",xlab="distance (m)")
legend(10,0.7,legend="1 m seed",lwd=2,lty=1,col='black',cex=0.7)
plot(x1,pd.20,xlim=c(0,100),type='l',lwd=2,ylab="probability density",xlab="distance (m)")
points(x1,pd.50,type='l',lwd=2,lty=2)
points(x1,pd.p,type='l',lwd=2,col='blue')
legend(60,6e-04,legend=c("20 m seed","50 m seed","pollen"),lwd=c(2,2,2),lty=c(1,2,1),col=c('black',"black","blue"),cex=0.7)
dev.off()

#probability of dispersing to each 20 m distance category (up to 2010 m)
div <- seq(101,20101,by=200)
Sd.prob.1 <- rep(0,100)
Sd.prob.20 <- rep(0,100)
Sd.prob.50 <- rep(0,100)
P.prob <- rep(0,100)
for(i in 1:100) {
  if(i==1){
  	  Sd.prob.1[i] <- sum(pd.1[1:div[i]])
  	  Sd.prob.20[i] <- sum(pd.20[1:div[i]])
  	  Sd.prob.50[i] <- sum(pd.50[1:div[i]])
  	  P.prob[i] <- sum(pd.p[1:div[i]])
  }
    if(i>1){
  	  Sd.prob.1[i] <- sum(pd.1[(div[i-1]+1):div[i]])
  	  Sd.prob.20[i] <- sum(pd.20[(div[i-1]+1):div[i]])
  	  Sd.prob.50[i] <- sum(pd.50[(div[i-1]+1):div[i]])
  	  P.prob[i] <- sum(pd.p[(div[i-1]+1):div[i]])
  }
	}
Sd.prob.1 <- Sd.prob.1/sum(Sd.prob.1)	
Sd.prob.20 <- Sd.prob.20/sum(Sd.prob.20)
Sd.prob.50 <- Sd.prob.50/sum(Sd.prob.50)
P.prob <- P.prob/sum(P.prob)
	
#number of cells per 20 m category
cell.d.num <- 1
for(i in 2:100){
	cell.d.num <- c(cell.d.num,8*(i-1))
}
#probability of dispersing to individual cells within 20 m distance categories (up to 2010 m)
Sd.prob.1.e <- Sd.prob.1/cell.d.num
Sd.prob.20.e <- Sd.prob.20/cell.d.num
Sd.prob.50.e <- Sd.prob.50/cell.d.num
P.prob.e <- P.prob/cell.d.num

#cell "coordinates"
X <- numeric(0)
Y <- numeric(0)
for(i in 1:100){
	X <- c(X,rep(i,50))
	Y <- c(Y,seq(1,50,by=1))
	}
cell.coord <- cbind(KK,X,Y)

#distance matrix
# to calculate distances
distmat <- function(x1,y1,x2,y2){
    xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)
    yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
    d <- t(sqrt(xd + yd)) 
    return(d)
}

cell.dist <- distmat(X,Y,X,Y)

##Seed dispersal matrix
SD.mat.20 <- matrix(0,length(KK),length(KK))
  for(i in KK){
  	  SD.mat.20[i,i]<-Sd.prob.20.e[1]
  	for(d in 2:100){
  		a <- which(cell.dist[i,]>=(d-1) & cell.dist[i,]<d)
  		SD.mat.20[i,a]<-Sd.prob.20.e[d]
  	}
  }
  
SD.mat.1 <- matrix(0,length(KK),length(KK))
  for(i in KK){
  	  SD.mat.1[i,i]<-Sd.prob.1.e[1]
  	for(d in 2:100){
  		a <- which(cell.dist[i,]>=(d-1) & cell.dist[i,]<d)
  		SD.mat.1[i,a]<-Sd.prob.1.e[d]
  	}
  }  
  
SD.mat.50 <- matrix(0,length(KK),length(KK))
  for(i in KK){
  	  SD.mat.50[i,i]<-Sd.prob.50.e[1]
  	for(d in 2:100){
  		a <- which(cell.dist[i,]>=(d-1) & cell.dist[i,]<d)
  		SD.mat.50[i,a]<-Sd.prob.50.e[d]
  	}
  }

##Pollen dispersal matrix  
 Pol.mat <- matrix(0,length(KK),length(KK))
  for(i in KK){
  	  Pol.mat[i,i]<P.prob.e[1]
  	for(d in 2:100){
  		a <- which(cell.dist[i,]>=(d-1) & cell.dist[i,]<d)
  		Pol.mat[i,a]<-P.prob.e[d]
  	}
  }

#jpeg(filename="SeedDisp_ann_20.jpg",width=9,height=7, units='in', res=500)
jpeg(filename="SeedDisp_ann_1.jpg",width=9,height=7, units='in', res=500)
#Seeds <- SD.mat.20[3025,]*5000
Seeds <- SD.mat.1[3025,]*5000

Seed.Mat <- matrix(Seeds,ncol=Plot.width,nrow=Plot.length,byrow=T) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=Seed.Mat,col=c("white","yellow","gold1","gold2","orange","dark orange","red"),breaks=c(0,0.00001,0.001,0.1,10,100,1000,5000),xlab="Distance",ylab="Distance",main="Seeds (out of 5000)") 
 legend(1500,950,legend=c("<0.00001","<0.001","<0.1","<10","<100","<1000","<5000"), fill=c("white","yellow","gold1","gold2","orange","dark orange","red") ,cex=0.7) 
dev.off() 

jpeg(filename="SeedDisp_Fig.jpg",width=9,height=7, units='in', res=500)
Seeds2 <- SD.mat.20[3025,]*5000
Seeds1 <- SD.mat.1[3025,]*5000

Seed.Mat1 <- matrix(Seeds1,ncol=Plot.width,nrow=Plot.length,byrow=T) 
Seed.Mat2 <- matrix(Seeds2,ncol=Plot.width,nrow=Plot.length,byrow=T)

par(mfrow=c(2,1), mar = c(4.5,4,2.5,1.5)) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=Seed.Mat1,col=c("white","yellow","gold1","gold2","orange","dark orange","red"),breaks=c(0,0.00001,0.001,0.1,10,100,1000,5000),xlab="Distance(m)",ylab="Distance (m)",main="1 m mean dispersal",cex.main = 1) 
 legend(1500,950,legend=c("<0.00001","<0.001","<0.1","<10","<100","<1000","<5000"), fill=c("white","yellow","gold1","gold2","orange","dark orange","red") ,cex=0.7) 
 
image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=Seed.Mat2,col=c("white","yellow","gold1","gold2","orange","dark orange","red"),breaks=c(0,0.00001,0.001,0.1,10,100,1000,5000),xlab="Distance(m)",ylab="Distance (m)",main="20 m mean dispersal",cex.main = 1) 
 
dev.off() 

save.image(file="AnnMod_LA_Setup.RData")