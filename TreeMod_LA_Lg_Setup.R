#50 x 200 20mx20m plot landscape w/10 climates.  
# Columns 1-20 climate 10, col 21-40 climate 9, col 41-60 climate 8; col 61-80 climate 7; col 81-100 climate 6; col 101-120 climate 5; col 121-140 climate 4; col 141-160 climate 3; col 161-180 climate 2; and col 181-200 climate 1.

#For genotype-specific simulation, varying proportions of each genotype by climate band.

#Output: Map of initial allele distribution, text file of population numbers, dispersal kernel graphs, dispersal matrices saved in R file 

##############Setup plot
Plot.width <- 50
Plot.length <- 200
C.width <- Plot.length/10
CellNum <- (Plot.width*Plot.length)

 #Habitat cells
   KK <- 1:CellNum
   #habitat cell quality - climate, given by 1 (cold) to 10 (hot)
kq <- c(rep(10,Plot.width*C.width),rep(9,Plot.width*C.width),rep(8,Plot.width*C.width),rep(7,Plot.width*C.width),rep(6,Plot.width*C.width),rep(5,Plot.width*C.width),rep(4,Plot.width*C.width),rep(3,Plot.width*C.width),rep(2,Plot.width*C.width),rep(1,Plot.width*C.width)) 

   
   #Cell indexing - form into rectangle
   index <- cbind(row=rep(1:Plot.length,each=Plot.width),column=rep(1:Plot.width,Plot.length)) #xy coord
   
 clim.state <- matrix(kq,ncol=Plot.width,nrow=Plot.length,byrow=T)
 
 
########Tree population
 #Size classes - 5 size classes
   JJ <- 1:5

   #Average basal area for each class
   BAJJ <- c(0.0001,0.003,0.07,0.14,0.5)

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
   	if(kq[k]==10)  N.i[k,2,] <- c(0,0,0,0,10,30)  #saplings -40
   	if(kq[k]==9)   N.i[k,2,] <- c(0,0,2,2,26,50)  #saplings - 80
   	if(kq[k]==8)   N.i[k,2,] <- c(0,10,15,15,60,100)  #saplings - 200
   	if(kq[k]==7)  N.i[k,2,] <- c(0,20,25,25,100,30)  #saplings - 200
     	if(kq[k]==6)   N.i[k,2,] <- c(10,40,45,45,50,10)  #saplings - 200
      	if(kq[k]==5)   N.i[k,2,] <- c(30,100,25,25,20,0)  #saplings - 200
     	if(kq[k]==4)   N.i[k,2,] <- c(100,60,15,15,10,0)  #saplings - 200
   	if(kq[k]==3)   N.i[k,2,] <- c(50,26,2,2,0,0)  #saplings - 80
   	if(kq[k]==2)   N.i[k,2,] <- c(30,10,0,0,0,0)  #saplings - 40
   }
   
    ##where are the genotypes (saplings)?  
   pop.state <- rep(0,length(KK))
   prop1 <- rep(0,length(KK)) #proportion of A1
   prop2 <- rep(0,length(KK)) #proportion of A2
   prop3 <- rep(0,length(KK)) #proportion of A3
   
   for(k in KK){
   	if(sum(N.i[k,2,II])>0){
   		 prop1[k] <-  (sum(0.5*N.i[k,2,2])+sum(0.5*N.i[k,2,4])+sum(N.i[k,2,1]))/sum(N.i[k,2,II]) 
   		 prop2[k] <-  (sum(0.5*N.i[k,2,5])+sum(0.5*N.i[k,2,2])+sum(N.i[k,2,3]))/sum(N.i[k,2,II]) 
   		 prop3[k] <-  (sum(0.5*N.i[k,2,4])+sum(0.5*N.i[k,2,5])+sum(N.i[k,2,6]))/sum(N.i[k,2,II]) 
   		  	
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
     
 jpeg(filename="Initial_Tree_LA_Lg_alleles.jpg",width=9,height=7, units='in', res=500)
par(mfrow=c(2,2),mar=c(4,4,3,2))  
 image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop1.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A1") 
 legend(50,950,legend=c("<0.1","0.1-0.3","0.3-0.5","0.5-0.7",">0.7"), fill=c("white","yellow","gold2","green2","green4") ,cex=0.7) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop2.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A2") 
 
image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=prop3.mat,col=c("white","yellow","gold2","green2","green4"),breaks=c(0,0.1,0.3,0.5,0.7,1),xlab="Distance",ylab="Distance",main="A3") 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=clim.state ,col=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5),xlab="Distance",ylab="Distance",main="Climate") 
legend(50,950,legend=c("1","2","3","4","5","6","7","8","9","10"), fill=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),cex=0.7)  
dev.off()    


####### Vegetative processes
Ideal <- c(4,5,6,6,7,8)   #Ideal climate for each genotype 

max.fec <- c(0,0,150,1000,5000)  #maximum seed production for each size class
max.pol <- c(0,0,800,10000,100000)  #maximum pollen production for each size class
max.germ <- 0.6				#maximum germination (for now no seed bank)
max.surv <- c(0.6,0.85,0.96,0.98,0.99) #maximum survival for each size class
max.trans <-c(0.15,0.11,0.06,0.02)  #maximum transition prob for each size class

#BA.eff <- c(0.06,0.03,0.015,0.01,0) #basal area effect on survival and transition 0.0001
BA.eff <- c(0.1,0.07,0.02,0.01,0) #basal area effect on survival and transition 0.0001

SFclim.1.6 <-c(0.9,0.5,0.3,0.15,0,0,0,0,0)  #survival/fecundity effect of climate, genotypes 1 & 6; ideal, 1-6 steps
SFclim.2.5 <- c(1,0.7,0.4,0.25,0.1,0,0,0,0) #survival/fecundity effect of climate, genotypes 2-5; ideal, 1-6 steps
#SFclim.1.6 <- c(0.9,0.75,0.55,0.4,0.3,0.2,0.1,0,0) #survival/fecundity effect of climate, genotypes 1 & 6; ideal, 1-6 steps
#SFclim.2.5 <- c(1,0.8,0.6,0.45,0.35,0.2,0.1,0,0) #survival/fecundity effect of climate, genotypes 2-5; ideal, 1-6 steps

TGclim <- c(1,0.7,0.5,0.2,0.1,0,0,0,0) #transition/germination effect of climate; ideal, 1 step-6 steps
#TGclim <- c(1,0.9,0.8,0.6,0.4,0.2,0.1,0,0) #transition/germination effect of climate; ideal, 1 step-6 steps


##########Dispersal
###2D-t seed dispersal
#mean dispersal = 40
u.40 <- (40*2/pi)^2 #648.5
#mean dispersal = 80
u.80 <- (80*2/pi)^2 #2593.8

#pollen mean dispersal = 180
u.p <- (180*2/pi)^2 #13131.2

x1 <- seq(0,2010,by=0.1)
pd.40 <- 1/(pi*u.40*(1+((x1^2)/u.40))^2)
pd.80 <- 1/(pi*u.80*(1+((x1^2)/u.80))^2)
pd.p <- 1/(pi*u.p*(1+((x1^2)/u.p))^2)

#probability of dispersing to each 20 m distance category (up to 2010 m)
div <- seq(101,20101,by=200)
Sd.prob.40 <- rep(0,100)
Sd.prob.80 <- rep(0,100)
P.prob <- rep(0,100)
for(i in 1:100) {
  if(i==1){
  	  Sd.prob.40[i] <- sum(pd.40[1:div[i]])
  	  Sd.prob.80[i] <- sum(pd.80[1:div[i]])
  	  P.prob[i] <- sum(pd.p[1:div[i]])
  }
    if(i>1){
  	  Sd.prob.40[i] <- sum(pd.40[(div[i-1]+1):div[i]])
  	  Sd.prob.80[i] <- sum(pd.80[(div[i-1]+1):div[i]])
  	  P.prob[i] <- sum(pd.p[(div[i-1]+1):div[i]])
  }
	}
Sd.prob.40 <- Sd.prob.40/sum(Sd.prob.40)	
Sd.prob.80 <- Sd.prob.80/sum(Sd.prob.80)
P.prob <- P.prob/sum(P.prob)
	
#number of cells per 20 m category
cell.d.num <- 1
for(i in 2:100){
	cell.d.num <- c(cell.d.num,8*(i-1))
}
#probability of dispersing to individual cells within 20 m distance categories (up to 2010 m)
Sd.prob.40.e <- Sd.prob.40/cell.d.num
Sd.prob.80.e <- Sd.prob.80/cell.d.num
P.prob.e <- P.prob/cell.d.num

#cell "coordinates"
X <- numeric(0)
Y <- numeric(0)
for(i in 1:200){
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
SD.mat.40 <- matrix(0,length(KK),length(KK))
  for(i in KK){
  	  SD.mat.40[i,i]<-Sd.prob.40.e[1]
  	for(d in 2:100){
  		a <- which(cell.dist[i,]>=(d-1) & cell.dist[i,]<d)
  		SD.mat.40[i,a]<-Sd.prob.40.e[d]
  	}
  }
  
  SD.mat.80 <- matrix(0,length(KK),length(KK))
  for(i in KK){
  	  SD.mat.80[i,i]<-Sd.prob.80.e[1]
  	for(d in 2:100){
  		a <- which(cell.dist[i,]>=(d-1) & cell.dist[i,]<d)
  		SD.mat.80[i,a]<-Sd.prob.80.e[d]
  	}
  }
  
 Pol.mat <- matrix(0,length(KK),length(KK))
  for(i in KK){
  	  Pol.mat[i,i]<P.prob.e[1]
  	for(d in 2:100){
  		a <- which(cell.dist[i,]>=(d-1) & cell.dist[i,]<d)
  		Pol.mat[i,a]<-P.prob.e[d]
  	}
  }


save.image(file="TreeMod_LA_Lg_Setup.RData")