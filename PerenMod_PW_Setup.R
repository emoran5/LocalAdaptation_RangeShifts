#50 x 100 20mx20m plot landscape w/10 climates.  
# Columns 1-10 climate 10, col 11-20 climate 9, col 21-30 climate 8; col 31-40 climate 7; col 41-50 climate 6; col 51-60 climate 5; col 61-70 climate 4; col 71-80 climate 3; col 81-90 climate 2; and col 91-100 climate 1.

#Output: Map of initial distribution, text file of population numbers, dispersal kernel graphs, dispersal matrices saved in R file 

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
 
 
########Population
 #Size classes - 4 size classes.  Seedling, small adult, medium adult, large adult
   JJ <- 1:4

   #Average size for each class (m2)
   AJJ <- c(0.004,0.06,0.25,2)

   #Genotypes - just one 
   II <- 1  
   
      #Initial population structure
   n.entry <- length(KK)*length(II)*length(JJ) 
   #habitat cells (5000) x genotypes (6) x size classes (5)
   N.i <- rep(0,n.entry)
   dim(N.i) <- c(length(KK),length(JJ),length(II)) #5000 x 5 matrix, set of 6
   
   ##where they start out
   for(k in KK){
   #Climate 1...10
   	if(kq[k]==10) N.i[k,1,] <- 400  #seedlings 
   	if(kq[k]==9)  N.i[k,1,] <- 800   
   	if(kq[k]==8)  N.i[k,1,] <- 2000  
   	if(kq[k]==7)  N.i[k,1,] <- 2000  
   	if(kq[k]==6)  N.i[k,1,] <- 2000  
	if(kq[k]==5)  N.i[k,1,] <- 2000   
   	if(kq[k]==4)  N.i[k,1,] <- 2000  
	if(kq[k]==3)  N.i[k,1,] <- 800 
	if(kq[k]==2)  N.i[k,1,] <- 400   
   }   
    ##where are the genotypes (seedlings)?  
   pop.state <- rep(0,length(KK))
   
   for(k in KK){
   	if(sum(N.i[k,1,II])>0) pop.state[k] <-  1   		  	
   }
   		 
   mat.pop.state <- matrix(pop.state,ncol=Plot.width,nrow=Plot.length,byrow=T) 
     
 jpeg(filename="Initial_Peren_P.jpg",width=9,height=7, units='in', res=500)
par(mfrow=c(2,2),mar=c(4,4,3,2))  
 image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=mat.pop.state,col=c("white","green4"),breaks=c(0,0.7,1),xlab="Distance",ylab="Distance",main="A1") 
 legend(50,950,legend=c("absent","present"), fill=c("white","green4") ,cex=0.7) 

image(x=c((1:Plot.length)*20),y=c((1:Plot.width)*20),z=clim.state ,col=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5),xlab="Distance",ylab="Distance",main="Climate") 
legend(50,950,legend=c("1","2","3","4","5","6","7","8","9","10"), fill=c("mediumseagreen","seagreen1","yellowgreen","yellow","gold1","gold3","orange","darkorange2","red1","red3"),cex=0.7)  
dev.off()    


#export as text file
Pop.i <- N.i
write.table(Pop.i,"InitialPop_Peren_P.txt") 


####### Vegetative processes
Ideal <- 6  #Ideal climate 

max.fec <- c(0,150,1000,5000)  #maximum seed production for each size class
max.pol <- c(0,800,10000,100000)  #maximum pollen production for each size class
max.germ <- 0.45				#maximum germination (for now no seed bank)
max.surv <- c(0.3,0.5,0.65,0.7) #maximum survival for each size class
max.trans <-c(0.3,0.15,0.1)  #maximum transition prob for each size class
 
CA.eff <- c(0.005,0.003,0.001,0.0005) #competitive area effect on survival and transition 

SFclim <- c(1,1,1,0.9,0.5,0.3,0.15,0,0) #survival/fecundity effect of climate, ideal, 1-8 steps

TGclim <- c(1,1,1,1,0.7,0.5,0.2,0.1,0) #transition/germination effect of climate; ideal, 1 step-8 steps


##########Dispersal
###2D-t seed dispersal
#mean dispersal = 20
u.20 <- (20*2/pi)^2 #162.1
#mean dispersal = 80
u.80 <- (80*2/pi)^2 #2593.8

#pollen mean dispersal = 100
u.p <- (100*2/pi)^2 #4052.847

x1 <- seq(0,2010,by=0.1)
pd.20 <- 1/(pi*u.20*(1+((x1^2)/u.20))^2)
pd.80 <- 1/(pi*u.80*(1+((x1^2)/u.80))^2)
pd.p <- 1/(pi*u.p*(1+((x1^2)/u.p))^2)

#probability of dispersing to each 20 m distance category (up to 2010 m)
div <- seq(101,20101,by=200)
Sd.prob.20 <- rep(0,100)
Sd.prob.80 <- rep(0,100)
P.prob <- rep(0,100)
for(i in 1:100) {
  if(i==1){
  	  Sd.prob.20[i] <- sum(pd.20[1:div[i]])
  	  Sd.prob.80[i] <- sum(pd.80[1:div[i]])
  	  P.prob[i] <- sum(pd.p[1:div[i]])
  }
    if(i>1){
  	  Sd.prob.20[i] <- sum(pd.20[(div[i-1]+1):div[i]])
  	  Sd.prob.80[i] <- sum(pd.80[(div[i-1]+1):div[i]])
  	  P.prob[i] <- sum(pd.p[(div[i-1]+1):div[i]])
  }
	}
Sd.prob.20 <- Sd.prob.20/sum(Sd.prob.20)	
Sd.prob.80 <- Sd.prob.80/sum(Sd.prob.80)
P.prob <- P.prob/sum(P.prob)
	
#number of cells per 20 m category
cell.d.num <- 1
for(i in 2:100){
	cell.d.num <- c(cell.d.num,8*(i-1))
}
#probability of dispersing to individual cells within 20 m distance categories (up to 2010 m)
Sd.prob.20.e <- Sd.prob.20/cell.d.num
Sd.prob.80.e <- Sd.prob.80/cell.d.num
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

save.image(file="PerenMod_PW_Setup.RData")