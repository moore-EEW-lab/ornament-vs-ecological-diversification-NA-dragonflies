# ten.se <- data.frame(log.size, log.size.se, log.rw, log.rw.se, log.color, log.col.se)
# ten.se $binom <- names(log.size)
# write.csv(ten.se, 'lib.spp.10se.csv', row.names = FALSE)

setwd('/Users/michaelmoore/desktop/Working Directory')

# packages 
library(phytools)
library(car)
library(MASS)
library(MuMIn)
library(geiger)


##  load in the data and phylogeny

ten.dat <- read.csv('lib.spp.10se.csv') # load data
rownames(ten.dat) <- ten.dat$binom 
phylo <- read.tree('lib.na.phylo.tre') 
dim(ten.dat) # check to make sure the data loaded in okay. Should output 72 7
summary(phylo) # check to make sure the phylogeny loaded in okay. Should be 72 tips and 71 nodes. max should = 120.317071 


############### patterns of diversification

###### patterns of size diversification
log.size <- ten.dat$log.size
log.size.se <- ten.dat$log.size.se
names(log.size) <- ten.dat$binom
names(log.size.se) <- ten.dat$binom

size.BM <- fitContinuous(phylo, log.size, SE = log.size.se, model = 'BM', control = list(hessian = TRUE)) # brownian motion
size.EB <- fitContinuous(phylo, log.size, SE = log.size.se, model = 'EB', bounds = list(a = c(min = -2, max = 0.5))) # early burst
size.white <- fitContinuous(phylo, log.size, SE = log.size.se, model = 'white') # white noise
size.delta <- fitContinuous(phylo, log.size, SE = log.size.se, model = 'delta') # delta model
size.OU <- fitContinuous(phylo, dat = log.size, log.size.se, model = 'OU') # OU model

aic.size <- c(size.BM$opt$aic, size.EB$opt$aic, size.white$opt$aic, size.delta$opt$aic, size.OU$opt$aic) # save the AIC values from the three models
names(aic.size) <- c('Brownian Motion', 'early burst', 'white', 'Delta', 'OU')

aic.size
# Brownian Motion     early burst           white           Delta              OU 
      # -46.14322       -44.30801        20.44525       -44.32610       -44.67232 

###### patterns of rw diversification
log.rw <- ten.dat$log.rw
log.rw.se <- ten.dat$log.rw.se
names(log.rw) <- ten.dat$binom
names(log.rw.se) <- ten.dat$binom

rw.BM <- fitContinuous(phylo, log.rw, SE = log.rw.se, model = 'BM') # brownian motion model
rw.EB <- fitContinuous(phylo, log.rw, SE = log.rw.se, model = 'EB', control = list(hessian = TRUE)) # early burst
rw.white <- fitContinuous(phylo, log.rw, SE = log.rw.se, model = 'white') # white
rw.delta <- fitContinuous(phylo, log.rw, SE = log.rw.se, model = 'delta') # delta model
rw.OU <- fitContinuous(phylo, log.rw, SE = log.rw.se, model = 'OU') #OU

aic.rw <- c(rw.BM$opt$aic, rw.EB$opt$aic, rw.white$opt$aic, rw.delta$opt$aic, rw.OU$opt$aic) 
names(aic.rw) <- c('Brownian Motion', 'early burst', 'white', 'Delta', 'OU') 

aic.rw
# Brownian Motion     early burst           white           Delta              OU 
      # -146.9629       -151.7089       -114.1602       -148.7965       -151.2378 

###### patterns of color diversification

log.color <- ten.dat$log.color
log.col.se <- ten.dat$log.col.se
names(log.color) <- ten.dat$binom
names(log.col.se) <- ten.dat$binom

color.BM <- fitContinuous(phylo, log.color, SE = log.col.se, model = 'BM', control = list(hessian = TRUE))
color.EB <- fitContinuous(phylo, log.color, SE = log.col.se, model = 'EB', bounds = list(a = c(min = -1, max = 0)))
color.white <-fitContinuous(phylo, log.color, SE = log.col.se, model = 'white')
color.delta <- fitContinuous(phylo, log.color, SE = log.col.se, model = 'delta', bounds = list(delta = c(min = exp(-500), max = 15))) 
color.OU <- fitContinuous(phy = phylo, dat = log.color, model = "OU", SE = log.col.se)

aic.color <- c(color.BM$opt$aic, color.EB$opt$aic, color.white$opt$aic, color.delta$opt$aic, color.OU$opt$aic) 
names(aic.color) <- c('Brownian Motion', 'early burst', 'white', 'Delta', 'OU')

aic.color
# Brownian Motion     early burst           white           Delta              OU 
       # 356.9669        358.9669        340.1961        335.2425        334.8998 

#### Next, compare rates of phenotypic evolution using the Adams (2013) likelihood ratio approach. Need to run his custom's scripts first. 
########## start here
CompareRates.multTrait <- function(phy, x, TraitCov=T, ms.err=NULL, ms.cov=NULL){
  #Compares LLik of R-matrix vs. LLik of R-matrix with constrained diagonal
  
  #TraitCov = TRUE assumes covariation among traits (default)
  #ms.err allows the incorporation of within-species measurement error. Input is a matrix of species (rows) by within-species variation for each trait (columns).
  #ms.cov allows the incorporation of within-species covariation between traits. Input is a matrix of species (rows) by within-species covariation for each pair of traits (columns). These must be provided in a specific order, beginning with covariation between trait 1 and the rest, then trait 2 and the rest, etc. For instance, for 4 traits, the columns are: cov_12, cov_13, cov_14, cov_23, cov_24 cov_34.
  
  #Some calculations adapted from 'evol.vcv' in phytools (Revell, 2012)
  
  library(MASS)
  x<-as.matrix(x)
  N<-nrow(x)
  p<-ncol(x)
  C<-vcv.phylo(phy)
  C<-C[rownames(x),rownames(x)]
  if (is.matrix(ms.err)){    
    ms.err<-as.matrix(ms.err[rownames(x),])}
  if (is.matrix(ms.cov)){    
    ms.cov<-as.matrix(ms.cov[rownames(x),])}
  
  #Cholesky decomposition function for diagonal-constrained VCV
  build.chol<-function(b){
    c.mat<-matrix(0,nrow=p,ncol=p)
    c.mat[lower.tri(c.mat)] <- b[-1]  
    c.mat[p,p]<-exp(b[1])
    c.mat[1,1]<-sqrt(sum((c.mat[p,])^2))
    if(p>2){
      for (i in 2:(p-1)){
        c.mat[i,i]<-ifelse( (c.mat[1,1]^2-sum((c.mat[i,])^2) )>0,
                            sqrt(c.mat[1,1]^2-sum((c.mat[i,])^2)), 0)
      }}
    return(c.mat) 
  }
  
  #Fit Rate matrix for all traits: follows code of L. Revell (evol.vcv)
  a.obs<-colSums(solve(C))%*%x/sum(solve(C))   
  D<-matrix(0,N*p,p)
  for(i in 1:(N*p)) for(j in 1:p) if((j-1)*N<i&&i<=j*N) D[i,j]=1.0
  y<-as.matrix(as.vector(x))
  one<-matrix(1,N,1)
  R.obs<-t(x-one%*%a.obs)%*%solve(C)%*%(x-one%*%a.obs)/N
  if (TraitCov==F)    #for TraitCov = F
  { R.obs<-diag(diag(R.obs),p)  }
  #Calculate observed likelihood with or without measurement error
  LLik.obs<-ifelse(is.matrix(ms.err)==TRUE, 
                   -t(y-D%*%t(a.obs))%*%ginv((kronecker(R.obs,C)+ diag(as.vector(ms.err))))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                     determinant((kronecker(R.obs,C)+ diag(as.vector(ms.err))))$modulus[1]/2 , 
                   -t(y-D%*%t(a.obs))%*%ginv(kronecker(R.obs,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                     determinant(kronecker(R.obs,C))$modulus[1]/2
  ) 
  
  #Fit common rate for all traits; search over parameter space   
  sigma.mn<-min(diag(R.obs))   #reasonable start value for diagonal
  
  #Within-species measurement error matrix
  if(is.matrix(ms.err)){m.e<-diag(as.vector(ms.err))}
  
  #Within-species measurement error and trait covariation matrix
  if (is.matrix(ms.err) && is.matrix(ms.cov)){	
    within.spp<-cbind(ms.err,ms.cov)
    rc.label<-NULL
    for (i in 1:p){ rc.label<-rbind(rc.label,c(i,i)) }
    for (i in 1:p){
      for (j in 2:p){ if (i!=j && i<j){rc.label<-rbind(rc.label,c(i,j))} }}
    m.e<-NULL
    for (i in 1:p){
      tmp<-NULL
      for (j in 1:p){
        for (k in 1:nrow(rc.label)){
          if(setequal(c(i,j),rc.label[k,])==T) {tmp<-cbind(tmp,diag(within.spp[,k]))}
        }
      }
      m.e<-rbind(m.e,tmp)
    }
  }
  
  #likelihood optimizer for no trait covariation
  lik.covF<-function(sigma){  
    R<-R.obs
    diag(R)<-sigma
    LLik<-ifelse(is.matrix(ms.err)==TRUE, 
                 -t(y-D%*%t(a.obs))%*%ginv((kronecker(R,C)+ m.e))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant((kronecker(R,C)+ m.e))$modulus[1]/2 , 
                 -t(y-D%*%t(a.obs))%*%ginv(kronecker(R,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant(kronecker(R,C))$modulus[1]/2
    ) 
    if (LLik == -Inf) { LLikk <- -1e+10  }
    return(-LLik)
  }
  
  #likelihood optimizer with trait covariation
  lik.covT<-function(sigma){  
    low.chol<-build.chol(sigma)
    R<-low.chol%*%t(low.chol)
    
    LLik<-ifelse(is.matrix(ms.err)==TRUE, 
                 -t(y-D%*%t(a.obs))%*%ginv((kronecker(R,C)+ m.e))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant((kronecker(R,C)+ m.e))$modulus[1]/2 , 
                 -t(y-D%*%t(a.obs))%*%ginv(kronecker(R,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant(kronecker(R,C))$modulus[1]/2
    ) 
    if (LLik == -Inf)  {LLikk <- -1e+10  }
    return(-LLik)
  }
  
  ##Optimize for no trait covariation
  if (TraitCov==F)    
  { model1<-optim(sigma.mn,fn=lik.covF,method="Nelder-Mead", control = list(maxit = 500000, reltol = 5e-6))}
  ##Optimize with trait covariation
  R.offd<-rep(0,(p*(p-1)/2))
  if (TraitCov==T)  
  {model1<-optim(par=c(sigma.mn,R.offd),fn=lik.covT, method="Nelder-Mead", control = list(maxit = 500000, reltol = 5e-6))}
  
  #### Assemble R.constrained
  if (TraitCov==F){R.constr<-diag(model1$par,p)}
  if (TraitCov==T){  
    chol.mat<-build.chol(model1$par)
    R.constr<-chol.mat%*%t(chol.mat)}
  
  if(model1$convergence==0)
    message<-"Optimization has converged."
  else
    message<-"Optim may not have converged.  Consider changing start value or lower/upper limits."
  LRT<- (-2*((-model1$value-LLik.obs)))
  LRT.prob<-pchisq(LRT, (p-1),lower.tail=FALSE) #df = Nvar-1
  AIC.obs<- -2*LLik.obs+2*p+2*p #(2p twice: 1x for rates, 1x for anc. states)
  AIC.common<- -2*(-model1$value)+2+2*p #(2*1: for 1 rate 2p for anc. states)
  return(list(Robs=R.obs, Rconstrained=R.constr,Lobs=LLik.obs,Lconstrained=(-model1$value),LRTest=LRT,Prob=LRT.prob,
              AICc.obs=AIC.obs,AICc.constrained=AIC.common,optimmessage=message))   
 }

Rate.Correlation <- function(Tree, Data, Trait1, Trait2) {
 DataM<-as.matrix(Data)

  Rates<-CompareRates.multTrait(Tree, DataM[,c(which(colnames(DataM)==Trait1), which(colnames(DataM)==Trait2))])
  Rates

}


Trait.Correlation <- function(Tree, Data, Trait1, Trait2) {
T1 <- Data[,Trait1]
names(T1) <- rownames(Data)
T1.pic <- pic(T1, Tree)

T2 <- Data[,Trait2]
names(T2) <- rownames(Data)
T2.pic<-pic(T2, Tree)

lm.pic <- lm(T2.pic ~ T1.pic -1) 

plot(T1.pic, T2.pic, pch=19)
abline(lm.pic, col="red")

#summary(lm.pic)
Pvalue.pic <- anova(lm.pic)$'Pr(>F)'[1]
Pvalue.pic
}
######### end here


### to run this analysis, need to create matrices that includes our traits and the standard errors of the traits

adult.columns <- data.frame(log.size, log.rw, log.color) # combine our three named vectors into a "data frame"
adult.traits <- as.matrix(adult.columns) # then turn that dataframe into a matrix

ses.frame <- data.frame(log.size.se, log.rw.se, log.col.se) # combine our three named vectors for the standard errors into a data frame
adult.ses <- as.matrix(ses.frame) # then turn that dataframe into a matrix


## compare the rates with likelihood ratio test (Adams 2013)
# rate.test <- CompareRates.multTrait(phylo, adult.traits, TraitCov = T, ms.err = adult.ses, ms.cov = NULL)
# rate.test

# $Robs
               # log.size        log.rw    log.color
# log.size   5.957540e-04 -6.129719e-05 -0.001393426
# log.rw    -6.129719e-05  1.532961e-04  0.001512303
# log.color -1.393426e-03  1.512303e-03  0.175182638

# $Rconstrained
         # [,1]         [,2]         [,3]
# [1,] 1.000307 0.000000e+00 0.000000e+00
# [2,] 0.000000 1.000307e+00 1.533196e-05
# [3,] 0.000000 1.533196e-05 1.000307e+00

# $Lobs
# [1] -141.9236

# $Lconstrained
# [1] -626.3062

# $LRTest
# [1] 968.7653

# $Prob
# [1] 4.318144e-211

# $AICc.obs
# [1] 295.8471

# $AICc.constrained
# [1] 1260.612

# $optimmessage
# [1] "Optimization has converged."

