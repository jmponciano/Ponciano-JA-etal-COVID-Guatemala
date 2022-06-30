my.pca <- function(Xnp){
  
  Names <- row.names(Xnp)	
  dimX <- dim(Xnp);
  n    <- dimX[1];
  p    <- dimX[2];
  S    <- var(Xnp);
  Dm12       <- matrix(0,nrow=p,ncol=p);
  diag(Dm12) <- 1/sqrt(diag(S));
  R          <- Dm12%*%S%*%Dm12;
  
  x.bar      <- apply(Xnp,2,mean);
  Znp        <- matrix(0,nrow=n,ncol=p);
  for(i in 1:n){
    Znp[i,] <- Dm12%*%(Xnp[i,]-x.bar);	
  }
  eigs.R <- eigen(R);
  eigs.S <- eigen(S);
  
  princomp.scoresR <- matrix(0,nrow=n,ncol=p); # Componentes principales ('scores')
  corr.pcompi.varkR <- matrix(0,nrow=p,ncol=p);
  prop.varsR        <- rep(0,p);
  
  princomp.scoresS <- matrix(0,nrow=n,ncol=p); # Componentes principales ('scores')
  corr.pcompi.varkS <- matrix(0,nrow=p,ncol=p);
  prop.varsS        <- rep(sum(eigs.S$values),p);
  diagS             <- diag(S);
  
  for(i in 1:p){
    
    princomp.scoresR[,i]  <- Znp%*%eigs.R$vectors[,i];
    corr.pcompi.varkR[,i] <- sqrt(eigs.R$values[i])*eigs.R$vectors[,i];
    prop.varsR[i]         <- eigs.R$values[i]/p;
    
    princomp.scoresS[,i]  <- Xnp%*%eigs.S$vectors[,i];
    corr.pcompi.varkS[,i] <- (sqrt(eigs.S$values[i])*eigs.S$vectors[,i])/sqrt(diagS);
    prop.varsS[i]         <- eigs.S$values[i]/prop.varsS[i];
    
  }
  
  # Row and col names of the resulting matrices:
  names.corrs <- list(colnames(Xnp),paste(rep("PrinComp",p), 1:p));
  names.scores <- list(Names,paste(rep("PrinComp",p), 1:p));
  
  dimnames(corr.pcompi.varkR) <- names.corrs;
  dimnames(corr.pcompi.varkS) <- names.corrs;
  
  dimnames(princomp.scoresR) <- names.scores;
  dimnames(princomp.scoresS) <- names.scores;            
  
  # The new coordinates of every obs. in each one of the principal component axis 
  #princomp.scoresR
  
  # The correlations of every variable with every principal component:  this is the stuff you want
  #corr.pcompi.varkR
  
  return(list(corr.pcomps.vars.R = corr.pcompi.varkR, corr.pcomps.vars.S= corr.pcompi.varkS, princomp.scoresR = princomp.scoresR,princomp.scoresS=princomp.scoresS, prop.varsR=prop.varsR, prop.varsS=prop.varsS))
  
}