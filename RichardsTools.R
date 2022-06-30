### helper functions for Richards fit analysis
### created:  6-08-21
### J.A.Ponciano


## richards likelihood function
richards.lhd <- function(guess, tseries, tvec){
	logl <- 0;
	r    <- 2./(1.+exp(guess[1]));
	a    <- 10/(1+exp(-guess[2])) +1;
	k    <- guess[3];
	tc   <- exp(guess[4]);
	no   <- tseries;
#	pred <- k/(1.+ exp(-r*(tvec-tc)))^(1/a);
	pred <- k/(1.+ a*exp(-r*(tvec-tc)))^(1/a);   #Wang et al,J.Theor.Bio (2012)
	ssq  <- sum((pred-tseries)^2);
	sig.hat <- sqrt(ssq/length(tseries));
	logl <- sum(dnorm(x=tseries, mean=pred, sd=sig.hat, log=TRUE))
#	if(is.finite(logl)!=TRUE){logl <- .Machine$double.xmax}
	nlog.like<- (-1)*logl
	return(nlog.like)	
}

## Richards growth fit predictions using estimated parameters  (par)

predictionW1<-function(par){
	pred.list <- list()
	tpred <- list()
	for(i in 1:7){
		tpred[[i]]  <- 1:(par$sup[i])
		k  <- par$k[i]
		a  <- par$a[i]
		r  <- par$r[i]
		tc <- par$tc[i]
#		pred.list[[i]]<- k/(1.+ exp(-r*(tpred[[i]]-tc)))^(1/a)
		pred.list[[i]]<- k/(1.+ a*exp(-r*(tpred[[i]]-tc)))^(1/a)		
		}
	return(c(tpred,pred.list))			
}

#### Richards fit function
richards.fit<-function(guess=c(1,2,22100,log(130)),tseries,deptnumber,pais="CR",figure=FALSE){
	if(pais=="CR"){
		deptos <-c("San Jose","Alajuela", "Cartago","Heredia", "Guanacaste","Puntarenas", "Limon")
	}else{deptos<-c("ALTA VERAPAZ", "BAJA VERAPAZ", "CHIMALTENANGO", "CHIQUIMULA", "EL PROGRESO",
          "ESCUINTLA", "GUATEMALA", "HUEHUETENANGO", "IZABAL", "JALAPA", "JUTIAPA", 
          "PETEN", "QUETZALTENANGO", "QUICHE", "RETALHULEU", "SACATEPEQUEZ", "SAN MARCOS",
          "SANTA ROSA", "SOLOLA", "SUCHITEPEQUEZ", "TOTONICAPAN", "ZACAPA")}
	
	i<-deptnumber
	depto.series  <- tseries[,1+i]
	tt <- 1:length(depto.series)
	ms.richards <- optim(par=guess,fn=richards.lhd, method="Nelder-Mead",tseries=depto.series,
	tvec= tt, hessian= TRUE)
	mse.list     <- ms.richards$par
	mval         <- ms.richards$val
#	hes          <- ms.richards$hessian 
#	fisher_info  <- solve(hes)
#	var.mls      <- diag(fisher_info)
#	std.mls      <- sqrt(var.mls)
#	z.alpha.half <- qnorm(p=0.975) 
#    upper        <- mse.list+z.alpha.half*std.mls
#    lower        <- mse.list-z.alpha.half*std.mls 
    convergencia <- ms.richards$convergence
    
	r.list       <- 2/(1+exp(mse.list[1]))
	a.list       <- 10/(1+exp(-mse.list[2])) +1
	k.list       <- mse.list[3]
	tc.list      <- exp(mse.list[4])
#	up.lim       <- c(2/(1+exp(upper[1])), 10/(1+exp(-upper[2])) +1, upper[3], exp(upper[4])) 
#	low.lim		 <- c(2/(1+exp(lower[1])), 10/(1+exp(-lower[2])) +1, lower[3], exp(lower[4]))
	params       <- c(mval,r.list, a.list, k.list, tc.list)
	names(params)<-c("likelihood","r","a","k","tc")
	
	names(k.list)  <- "k"
	names(a.list)  <- "a"
	names(r.list)  <- "r" 
	names(tc.list) <- "tc"

	if(figure==TRUE){
		
		tt <- 1:length(tseries[,1+i])
		k  <- k.list
		a  <- a.list
		r  <- r.list
		tc <- tc.list
#		pred.list<- k/(1.+ exp(-r*(tt-tc)))^(1/a)
		pred.list<- k/(1.+a*exp(-r*(tt-tc)))^(1/a)  # Wang et al

		ndays <- length(tseries[,1])
		tt    <- seq(1,ndays)
		
		jj<-deptnumber
		plot(tt,tseries[,jj+1],cex.lab=1.5,cex.axis=1.25,type="l",lty=1,lwd=2,bty="l",xlab="days",ylab="Cumulative infected",xlim=c(tt[1],tt[ndays])) 
		points(tt,pred.list, pch=16,cex=0.5) 
		mtext(deptos[jj], side=3,adj=0.1, font=3) 
	}
	 
	return(c(mval,r.list,a.list,k.list,tc.list,convergencia))
}



### cumulative time series for Costa Rica
### arguments are day0 (first day for fit analysis), file (csv data file with incidence cases) and ts ("cum" or "daily")

cum.CR<-function(day0,file,ts="cum"){
	deptos <-c("San Jose","Alajuela", "Cartago","Heredia", "Guanacaste","Puntarenas", "Limon")
	tseries<-read.csv(file)
	tseries$provincia<-gsub("San Jos\x8e","San Jose",tseries$provincia)
	tseries$provincia<-gsub("Lim\x97n","Limon",tseries$provincia)
	all.dates<-names(tseries[,-(1:4)])
	tseries <- tseries[,-(3:4)]
	tseries <- tseries[,-(1)]

	all.dates	 <-	 gsub("X","",all.dates)
	all.dates    <-  gsub("2020","20",all.dates)
	all.dates    <-  gsub("2021","21",all.dates)
	all.dates	 <-	 as.Date(all.dates,format="%d.%m.%y")
	len.dates    <-  length(all.dates)

	cum.ts       <-  data.frame(all.dates)
	ts.mat       <-  matrix(0,ncol=length(deptos),nrow=len.dates)
	daily.mat	 <-  matrix(0,ncol=length(deptos),nrow=len.dates)
	len.dep      <- length(deptos)

	for(i in  1:length(deptos)){
		one.depto      <-  tseries[tseries$provincia==deptos[i],]
		total.p.day    <-  mapply(sum,one.depto[,-(1:1)])
		ts.mat[,i]     <-  total.p.day
		}

	cum.ts <-  data.frame(all.dates,ts.mat)
	names(cum.ts)<-c("fechas",deptos)


	cum.ts    <- cum.ts[day0:len.dates,]   # starting from day0 in timeseries
	
	daily.mat[1,] <- ts.mat[1,]
	daily.mat[2:len.dates,] <- ts.mat[2:len.dates,]-ts.mat[1:(len.dates-1),]
	daily.ts <- data.frame(all.dates,daily.mat)
	names(daily.ts)<-c("fechas",deptos)
	
	daily.ts<-daily.ts[day0:len.dates,]
	if(ts=="cum"){final.ts<-cum.ts}else{final.ts<-daily.ts}
	
	return(final.ts)	
}

### cumulative time series for Guatemala
### arguments are day0 (first day for fit analysis) and file (csv data file with incidence cases)

cum.GT<-function(day0,file,ts="cum"){
	deptos<-c("ALTA VERAPAZ", "BAJA VERAPAZ", "CHIMALTENANGO", "CHIQUIMULA", "EL PROGRESO",
          "ESCUINTLA", "GUATEMALA", "HUEHUETENANGO", "IZABAL", "JALAPA", "JUTIAPA", 
          "PETEN", "QUETZALTENANGO", "QUICHE", "RETALHULEU", "SACATEPEQUEZ", "SAN MARCOS",
          "SANTA ROSA", "SOLOLA", "SUCHITEPEQUEZ", "TOTONICAPAN", "ZACAPA")
    
    tseries.all <- read.csv(file)
    #just deptos
    tseries      <-  tseries.all[,-(2:5)]
	all.dates	 <-  names(tseries)[-(1)]
	all.dates	 <-	 gsub("X","",all.dates)
	all.dates    <-  gsub("2020","20",all.dates)
	all.dates    <-  gsub("2021","21",all.dates)
	all.dates	 <-	 as.Date(all.dates,format="%y.%m.%d")
#	all.dates	 <-	 as.Date(all.dates,format="%m.%d.%y")
	len.dates    <-  length(all.dates)

	daily.ts     <-  data.frame(all.dates)
	ts.mat       <-  matrix(0,ncol=length(deptos),nrow=len.dates)

	for(i in  1:length(deptos)){
		one.depto      <-  tseries[tseries$departamento==deptos[i],]
		total.p.day    <-  mapply(sum,one.depto[,-(1:1)])
		ts.mat[,i]     <-  total.p.day
	}
	
	daily.ts <-  data.frame(all.dates,ts.mat)
	names(daily.ts)<-c("fechas",deptos)

	cum.ts  <- cumsum(daily.ts[,2:(length(deptos)+1)])
	cum.ts  <- data.frame(all.dates,cum.ts)
	names(cum.ts) <-c("fechas",deptos)
	
	cum.ts    <- cum.ts[day0:len.dates,]
	daily.ts  <- daily.ts[day0:len.dates,]
	if(ts=="cum"){final.ts<-cum.ts}else{final.ts<-daily.ts}
	
	return(final.ts)
}
