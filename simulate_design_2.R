simulateData <- function(n, ncp, p, sigma, sigmaNoise=0.1, ConcVarFact=0.8, n_min_peaks=3){
  meanPeakSigma <- sigma
  sigPeakSigma <- sigma / 4
  axis <- 1:p
  S <- matrix(0,ncp,p)
  C <- matrix(0,n,ncp)
  for(i in 1:ncp){
    npeaks <- 3+ceiling(10*runif(1))
    peakheights <- runif(npeaks)
    sigmas <- runif(npeaks) * sigPeakSigma + meanPeakSigma
    position <- runif(npeaks) * p
    for(j in 1:npeaks){
      S[i,] <- S[i,] + peakheights[j] * exp(-0.5 * ((axis - position[j]) / sigmas[j])^2)
    }
  }
  meanC <- sort(10^runif(ncp),decreasing = T)
  varC <- ConcVarFact * meanC * runif(ncp)
  for(i in 1:ncp){
    C[,i] <- rnorm(n = n,mean = meanC[i], sd = varC[i]/2)
  }
  X <- C%*%S;X <- X/max(abs(X))
  E <- matrix(rnorm(n*p,sd = sigmaNoise),nrow = n,ncol = p)
  X <- X*sqrt(1-sigmaNoise^2) + E
  list(X=X, C=C, S=S, E=E)
}

## Function to simulate data from *Design 2*
## The seed is chosen 1:100
simulateMulti <- function(seed=1,n=50,q=5,p1=500,p2=5000,
                          sigma1=0.05,sigma2=0.05,sigmaY=0.2,
                          ncpX=10,ncpXCom=5,ncpXYCom=3,plot=F){
  set.seed(seed)
  # ncpX for each X separately
  Data_1 <- simulateData(n=n, ncp=ncpX, p=p1, sigma=20, sigmaNoise = sigma1, ConcVarFact=0.8, n_min_peaks=5)
  Data_2 <- simulateData(n=n, ncp=ncpX, p=p2, sigma=30, sigmaNoise = sigma2, ConcVarFact=0.8, n_min_peaks=5)
  S1 <- Data_1$S;C1 <- Data_1$C;X1 <- Data_1$X
  S2 <- Data_2$S;C2 <- Data_2$C
  # ncpXCom in common
  C2[,1:ncpXCom] <- C1[,1:ncpXCom,drop=F]
  X2 <- C2%*%S2;X2 <- X2/max(abs(X2))
  E2 <- matrix(rnorm(n*p2,sd = sigma2),nrow = n,ncol = p2)
  X2 <- X2*sqrt(1-sigma2^2) + E2
  # Build Y on ncpXYCom components
  Y <- scale(Data_1$C[,1:ncpXYCom,drop=F])
  # Add extra variables and noise
  E_y <- matrix(rnorm(n*q,sd = sigmaY),nrow = n,ncol = q)
  if(q>ncpXYCom){
    Y <- cbind(Y,matrix(rnorm(n*(q-ncpXYCom)),nrow = n))*sqrt(1-sigmaY^2) + E_y
  }else{
    Y <- Y*sqrt(1-sigmaY^2) + E_y
  }
  if(plot){
    layout(matrix(c(1,2,2,3,3,3),nrow = 2,byrow = T))
    matplot(t(X1),lty=1,type="l")
    matplot(t(X2),lty=1,type="l")
    corXY <- cor(cbind(X1,X2),Y)
    matplot(corXY,type="l")
  }
  list(Xs=list(X1=X1,X2=X2),Y=Y,S=list(S1=S1,S2=S2,SXY=list(S1=S1[1:ncpXYCom,],S2=S2[1:ncpXYCom,])))
}


test <- function(){
  for(seed in 1:100){
    cat("\n_________________________________________________\n")
    cat(seed)
    datas <- simulateMulti(seed=seed,sigma1 = 0.05,sigma2 = 0.05,
                           sigmaY = sigma_y,plot=F)
    Xs <- datas$Xs;X_c <- cbind(Xs$X1,Xs$X2);Y <- datas$Y
    S_common <- cbind(datas$S$SXY$S1,datas$S$SXY$S2)
    ### Start the different methodologies on the chosen data set
  }

}
