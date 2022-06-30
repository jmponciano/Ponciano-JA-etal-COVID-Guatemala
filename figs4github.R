source("RichardsTools.R")
library(stringr)

## loading time series and demographic parameters
dia0 <- 50
cum.ts    <-cum.GT(day0=dia0,file="MediaMovilGT.csv")  # cumulative time series starting from day 50 
daily.ts  <-cum.GT(day0=dia0,file="MediaMovilGT.csv","daily") # daily time series starting from day 50

deptos    <- str_to_title(names(cum.ts)[2:23])  # list of provincias names
poblacion <- c(1215038,299476,615776,415063,176632,733181,3015081,1170669,408688,342923,488395,
               545600,799101,949262,326828,330469,1032277,396607,421583,554695,418569,245374)
trabajadores <-c(336847,93172,204163,107353,55771,212929,1239858,288481,123879,110893,
                 145024,152071,249657,246353,92157,132056,217281,124941,131977,165189,129263,73927)

estudiantes <-c(316747,79730,156220,100329,47074,185476,884829,256498,105971,90014,131359,147785,212824,
                225650,89118,89281,264774,108642,107405,144726,96538,60418)

cortesOlas  <- read.csv(file="cortesOlasGT.csv")

indicadores <- read.csv("IndicadoresGT2.csv")


##### fit analysis for First wave

pfinWave<- cortesOlas$F1-dia0  # data from media movil series starts at day 14
pinfl   <- cortesOlas$P1-dia0

par.mat<-matrix(0,nrow=22,ncol=4)
std.mat<-matrix(0,nrow=22,ncol=4)

for(n in 1:2){
  
  prov  <- n    ### departamento number (1 to 22)
  
  tf       <- pfinWave[prov]   # tf is the duration in days of the first wave. 
  cum.ts1  <- cum.ts[1:tf,]    # truncating the time series
  daily.ts1<- daily.ts[1:tf,]
  beta  <- 0.06                # inverse of temperature (simmulated annealing parameter)
  carr.cap <- tail(cum.ts1,1)[prov+1]
  
  result   <- richards.fit(guess=c(2,2,carr.cap,log(pinfl[prov])),tseries=cum.ts1,deptnumber=prov,pais="GT",figure=FALSE)
  
  
  ### first metropolis
  current_mval <- result[1]
  r.vec  <- result[2]
  a.vec  <- result[3]
  tc.vec <- tc <- result[5]
  k.vec  <- carr.cap <- result[4]
  lk.means <-c()
  
  beta.list<-seq(0.06,0.36,0.01)  #beta list for the annealing schedule
  print(result[1])
  for(beta in beta.list){
    lkhd.vec <-c()
    for(k in 1:300){
      delta <-carr.cap*runif(1,-0.05,0.05)
      result <-richards.fit(c(2,2,carr.cap+delta,log(tc)),tseries=cum.ts1,deptnumber=prov,pais="GT",figure=FALSE)
      delta_mval <-result[1]
      r <- runif(1,0,1)
      if( (r< exp(-beta*(delta_mval-current_mval))) ){
        current_mval<-delta_mval
        carr.cap<-carr.cap+delta
        result.new<-richards.fit(c(2,2,carr.cap,log(tc)),tseries=cum.ts1,deptnumber=prov,pais="GT",figure=FALSE)
        #print(result.new[6])
        r.vec <-append(r.vec,result.new[2])
        a.vec <-append(a.vec,result.new[3])
        k.vec <-append(k.vec,result.new[4])
        tc.vec<-append(tc.vec,result.new[5])
        lkhd.vec <- append(lkhd.vec,result.new[1])
      }
    }
    lk.means<-append(lk.means,mean(lkhd.vec))
  }
  k  <- mean(k.vec)
  a  <- mean(a.vec)
  r  <- mean(r.vec)
  tc <- mean(tc.vec)
  par.mat[prov,]<-c(a,r,k,tc)
  std.mat[prov,]<-c(sd(a.vec),sd(r.vec),sd(k.vec),sd(tc.vec))
  
  ### Plot with mean values of parameters
  tdays  <- length(cum.ts[,1])
  w1days <- length(cum.ts1[,1])
  
  td    <- seq(1,tdays)
  tt    <- seq(1,w1days)
  
  pred.list <- k/(1.+ a*exp(-r*(tt-tc)))^(1/a)
  wave      <- k*r*(1+a*exp(-r*(tt-tc)))^(-1/a)/(a+ exp(r*(tt-tc)))
  
  jj<-prov
  plot(td,daily.ts[,jj+1],cex.lab=1.5,cex.axis=1.25,type="l",lty=2,lwd=1,bty="l",xlab="days",ylab="Cumulative infected",xlim=c(td[1],td[tdays])) 
  points(tt,wave,pch=19,cex=0.5)
  mtext(deptos[jj], side=3,adj=0.1, font=3) 
}



####### principal component analysis


source("mypca.R")
library(MASS)
params1  <- read.csv(file = "W1paramsGT.csv")

mat  <- as.matrix(params1[,2:5])
deptos <- str_to_title(params1[,1])
rownames(mat) <- deptos
pca.out <- my.pca(mat)

# Proporcion de la varianza explicada (acumulada)
cum.var.prop <- cumsum(pca.out$prop.varsR)
# [1] 0.4816064 0.7672603 0.9690158 1.0000000

# 1st principal component scores
pca.out$princomp.scoresR[order(pca.out$princomp.scoresR[,1]),1]

p <- ncol(mat)

l1 <- matrix(pca.out$corr.pcomps.vars.R[,1],nrow=p,ncol=1) 
l2 <- matrix(pca.out$corr.pcomps.vars.R[,2],nrow=p,ncol=1) 


## mobility matrix (workers and students)
movilidad <- read.csv("MovilidadTrab-Est.csv")
loc.names <-movilidad[,1]
mov   <- movilidad[,-1]
nodes <- 22
adj   <- matrix(0,nrow=nodes,ncol=nodes,dimnames=list(loc.names,loc.names))
for(j in 1:nodes){
  adj[,j]      <-  mov[,j]
}


flow.pop<-matrix(0,nrow=nodes,ncol=nodes,dimnames=list(loc.names,loc.names))
for(jj in 1:nodes){
  flow.pop[jj,]<-adj[jj,]
}

## defining importation risk as No. of persons received by a locality divided by No. of persons leaving a locality of high risk
diag(flow.pop) <- rep(0,22)
leave.loc   <- sum(flow.pop[7,-7])+sum(flow.pop[17,-17])+sum(flow.pop[11,-11]) 
import.loc  <- colSums(flow.pop)
import.risk <- import.loc/leave.loc

#ordisurf
library(vegan)
library(MuMIn)

names(indicadores)
indicadores$ImpRisk <- 1/import.risk   # adding inverse of importation risk
pc.scores <- pca.out$princomp.scoresR[,1:2]

pca.df <- data.frame(PCI=-pc.scores[,1],PCII=-pc.scores[,2],indicadores)
rpca <- princomp(mat,cor=TRUE)
trial.surf4 <- ordisurf(rpca~ImpRisk,data=pca.df, plot=FALSE)




#########  Paper Fig. 1

params1  <- read.csv(file="W1paramsGT.csv")
tdays  <- length(cum.ts[,1])
td    <- seq(1,tdays)

par(mfrow=c(2,3),oma=c(1.5,1.5,2,1),mar=c(4,4,1,1))
for(jj in 13:18){
  a   <- params1[jj,2]
  r	<- params1[jj,3]
  k	<- params1[jj,4]
  tc	<- params1[jj,5]
  wave <-k*r*(1+a*exp(-r*(td-tc)))^(-1/a)/(a+ exp(r*(td-tc)))
  plot(td,daily.ts[,jj+1], cex.lab=1.,cex.axis=.9,type="l",lty=2,lwd=1,bty="l",xlab="",ylab="",col="gray",xlim=c(td[1],td[tdays])) 
  #	title(ylab="Daily infected",line=2)
  #	title(xlab="days",line=2)
  points(td,wave,pch=19,cex=0.3,col="blue")
  mtext(deptos[jj], side=3,adj=0.1, font=3,cex=0.6)	
}
mtext("Days",side=1,line=0,outer=TRUE,cex=1.)
mtext("Daily infected cases",side=2,line=0,outer=TRUE,cex=1.,las=0)

#### legend figure
par(fig=c(0,1,0,1), mar=c(0,0,0,0), oma=c(0,0,0,0), new=TRUE)
plot(0,0,type="n",bty="n",xaxt="n", yaxt="n")
legend("topright", legend=c("Richards Model", "Observed cases"), 
       lty=c(1,2), lwd=2, col=c("blue","grey"), bty="n",ncol=2)



### Paper Fig. 2

pos.vec <- rep(4,length(deptos))
pos.vec[2] <-3
pos.vec[4] <-3
pos.vec[5] <-2
pos.vec[8] <-2
pos.vec[10] <-2
pos.vec[12] <-1
pos.vec[13] <-1
pos.vec[15] <-3
pos.vec[20] <-1


par(oma=c(1,1,1,1), mar=c(4,4,1,1))
plot(-pca.out$princomp.scoresR[,1],-pca.out$princomp.scoresR[,2],pch=16, bty="n",
     xlab= "PC I", ylab= "PC II", xlim=c(-3,3.5), ylim=c(-3.5,3));
abline(h=0, lty=2)
abline(v=0, lty=2)
text(-pca.out$princomp.scoresR[,1],-pca.out$princomp.scoresR[,2],labels=deptos, 
     cex=0.7, pos=pos.vec)
arrows(x0=rep(0,p),y0=rep(0,p),x1=-l1,y1=-l2, length=0.08, lwd=2, col="blue");
text(-l1,-l2,labels=colnames(mat), pos=c(2,3,2,3),cex=0.85, col="blue")

plot(trial.surf4, add=TRUE, col="red",nlevels=16,labcex=1.05)
legend("topleft", legend=expression(paste("1/",rho)), lty=1,lwd=1,col="red", bty="n")


######### Paper Fig. 3
library(RColorBrewer)
count.cols <- c("blue","red","black","gray","cyan")

groups<-list(c(5,17,6,7),c(18,11,4),c(3,16,20,13),c(9,15,1,19,22),c(21,10,12),c(2,14,8)) ##obtained from PCA

par(mfrow=c(2,3),oma=c(1.5,1.5,2,1),mar=c(4,4,1,1))
for(iter in 1:length(groups)){
  plot(0,1,cex=0.,col="blue",xlab="",ylab="",xlim=c(0,405),ylim=c(0,1))
  ij<-1
  for(jj in groups[[iter]]){
    a    <- params1[jj,2]
    r	   <- params1[jj,3]
    k	   <- params1[jj,4]
    tc   <- params1[jj,5]
    wave <- k*r*(1+a*exp(-r*(td-tc)))^(-1/a)/(a+ exp(r*(td-tc)))
    max.wave <- max(wave)
    points(td,wave/max(wave),pch=19,cex=0.3,col=count.cols[ij])
    legend("topright", legend=deptos[groups[[iter]]], col=count.cols,pch=rep(19,13), cex=.80, bty="n")
    ij <- ij+1
  }
  mtext(paste0("G",iter),side=3,line=0.5, adj=0.5,font=1)
}
mtext("Days",side=1,line=0,outer=TRUE,cex=1.)
mtext("Daily infected cases",side=2,line=0,outer=TRUE,cex=1.,las=0)


###### Paper Fig. 5a
pos.vec[4] <-4
pos.vec[2] <-2
pos.vec[6] <-3
pos.vec[9] <-1
pos.vec[17]<-2
pos.vec[11]<-3
pos.vec[12]<-4
pos.vec[22]<-3
plot(import.risk,-pca.out$princomp.scoresR[,2],
     xlab= expression(rho), ylab="PC II" ,pch=16,cex=1.2, xlim=c(-0.05,0.8), ylim=c(-1,3))
abline(h=0, lty=2)
abline(v=0, lty=2)
text(import.risk, -pca.out$princomp.scoresR[,2],labels=deptos, 
     cex=0.7, pos=pos.vec)


##### Paper Fig. 5b
pos.vec[1] <-2
pos.vec[3] <-3
pos.vec[8] <-4
pos.vec[14] <-2
pos.vec[16] <-4
pos.vec[17] <-4
pos.vec[18] <-1
plot(1/import.risk,-pca.out$princomp.scoresR[,2],
     xlab= expression(paste("1/",rho)), ylab="PC II" ,pch=16,cex=1.2, xlim=c(0,30), ylim=c(-1,3))
abline(h=0, lty=2)
abline(v=0, lty=2)
text(1/import.risk, -pca.out$princomp.scoresR[,2],labels=deptos, 
     cex=0.7, pos=pos.vec)



##### Paper Fig. 6

library(cluster)
library(gam)
# without epicenters
pc1.short2 <- -pc.scores[-c(7,11,17,19,5),1]
pc2.short2 <- -pc.scores[-c(7,11,17,19,5),2]
ind.short2 <- indicadores[-c(7,11,17,19,5),]
df.short2 <- data.frame(PCI=pc1.short2,PCII=pc2.short2,ind.short2)
pc.scores.reduced   <- pc.scores[-c(7,11,17,19,5),]
import.risk.reduced <- import.risk[-c(7,11,17,19,5)]
deptos.reduced      <- deptos[-c(7,11,17,19,5)]
pos.vec.reduced     <- pos.vec[-c(7,11,17,19,5)]
colores <-brewer.pal(n=6,name="Paired")


GAM.M1 <- gam(PCII ~ ns(ImpRisk,2),data=df.short2)
gam.model <-GAM.M1

groups<-list(c(5),c(14,4),c(3,13,15,10),c(7,12,1,17),c(16,9,8),c(2,11,6)) ## matches fig wave classification for reduced list
puntos <-c(15,16,17,18,22,23)

pos.vec.reduced <- rep(4,length(deptos.reduced))
pos.vec.reduced[2] <-2
pos.vec.reduced[3] <-3
pos.vec.reduced[5] <-3
pos.vec.reduced[7] <-1
pos.vec.reduced[9] <-3
pos.vec.reduced[10] <-1
pos.vec.reduced[11] <-2
pos.vec.reduced[15] <-1
pos.vec.reduced[17] <-3

plot(gam.model, se=TRUE,
     xlab= expression(paste("1/",rho)), ylab="PC II" ,cex=1.2, xlim=c(0,30), ylim=c(-1,2),col="red")
abline(h=0, lty=2)
abline(v=0, lty=2)
text(1/import.risk.reduced, -pc.scores.reduced[,2],labels=deptos.reduced, 
     cex=0.7, pos=pos.vec.reduced)

for(i in 1:length(groups)){
  for(jk in groups[[i]]){
    points(1/import.risk.reduced[jk],-pc.scores.reduced[jk,2],pch=puntos[i],cex=1.2,col=colores[i])
  }
}
legend("topright",legend=c("G1","G2","G3","G4","G5","G6"), col=colores,pch=puntos,cex=0.9, ncol=1,bty="n")




