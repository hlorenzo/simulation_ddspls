get_design_1 <- function(n=50,sqrt_1_minus_sig2=0.99,p=1000,q=3){ # Structure
  alpha3 <- 1/sqrt(3)
  alpha2 <- 1/sqrt(2)
  repX <- 50
  A1 <- c(rep(alpha3,repX),rep(0,p-repX))
  A2 <- c(rep(0,repX),rep(alpha2,repX),rep(0,p-2*repX))
  A <- matrix(c(rep(A1,3),rep(A2,2)),nrow = 5,byrow = T)*sqrt_1_minus_sig2
  D1 <- c(rep(alpha3,1),rep(0,q-1))
  D2 <- c(rep(0,1),rep(alpha2,1),rep(0,q-2))
  D <- matrix(c(rep(D1,3),rep(D2,2)),nrow = 5,byrow = T)*sqrt_1_minus_sig2 # Observations
  d <- ncol(A)+nrow(A)+ncol(D)
  psi <- MASS::mvrnorm(n = n,mu = rep(0,d),Sigma = diag(d))
  phi <- psi[,1:nrow(A)]
  epsilonX_info <- psi[,nrow(A)+1:(2*repX)]*sqrt(1-sqrt_1_minus_sig2^2)
  epsilonX_noise <- psi[,nrow(A)+(2*repX)+1:(ncol(A)-2*repX)]
  epsilonY_info <- psi[,nrow(A)+ncol(A)+1:2]*sqrt(1-sqrt_1_minus_sig2^2)
  epsilonY_noise <- psi[,d]
  # X and Y
  X <- phi%*%A + cbind(epsilonX_info,epsilonX_noise)
  Y <- phi%*%D + cbind(epsilonY_info,epsilonY_noise)
  list(X=X,Y=Y)
}


nbs <- c(200,200,100,50,50)
ns <- c(25,50,100,200,400)
ids <- 1:500
NCORES <- 20
idd <- unlist(lapply(1:NCORES,function(ii){rep(ii,25)}))
`%my_do%` <- `%dopar%`
cl <- makeCluster(NCORES)
registerDoParallel(cl)
res <- foreach(ii=1:NCORES,
               .packages = "ddsPLS2",
               .combine='c',
               .multicombine=TRUE) %my_do% {
                 models_lambda0 <- list()
                 iiii <- c(1:500)[which(idd==ii)]
                 for(i in iiii){
                   x <- Xs[[i]]$x
                   y <- Ys[[i]]
                   n <- nrow(y)
                   models_lambda0[[i]] <- ddsPLS(x,y,verbose = F,n_B = nbs[which(ns==n)])
                 }
                 models_lambda0
               }
stopCluster(cl)

res2 <- res[which(unlist(lapply(res,length))>0)]

set.seed(1)
data_test <- get_design_1(n=1000)
MstdMSE <- matrix(NA,500,3)
for(i in 1:500){
  y_est <- predict(res2[[i]],data_test$X)
  MstdMSE[i,] <- colMeans((y_est$y_est-data_test$Y)^2)/apply(data_test$Y,2,sd)
}

ns <- unlist(lapply(res2,function(r){nrow(r$X)}))
Rs <- unlist(lapply(res2,function(r){r$R}))
Q2s <- unlist(lapply(res2,function(r){r$Q2[r$R] }))
selX <- lapply(res2,function(r){r$Selection$X })
selY <- lapply(res2,function(r){r$Selection$Y })
nx <- unlist(lapply(selX,length))
tprx <- unlist(lapply(selX,function(ii){length(which(ii %in% 1:100))/100}))
fprx <- unlist(lapply(selX,function(ii){length(!which(ii %in% 1:100))/900}))
ny <- unlist(lapply(selY,length))
dd <- data.frame(list(n=ns,R=Rs,Q2=Q2s,nx=nx,ny=ny,tprx=tprx,fprx=fprx,MstdMSE=MstdMSE))
par(mfrow=c(1,3));boxplot(MstdMSE.1~n,dd,ylim=c(0,1.4));abline(h=0,lty=3);boxplot(MstdMSE.2~n,dd,ylim=c(0,1.4));abline(h=0,lty=3);boxplot(MstdMSE.3~n,dd,ylim=c(0,1.4));abline(h=1,lty=3);
boxplot(Q2~n,dd);abline(h=0.6534)
boxplot(nx~n,dd);abline(h=100)
boxplot(ny~n,dd);abline(h=2)
boxplot(tprx~n,dd)
boxplot(fprx~n,dd)

for(i in 1:7){
  models <- ALL_FUCKING_MODELS[names(ALL_FUCKING_MODELS)[i]]
  save(models,file = paste("../data/SADM_simu_5_1_ALL_PREVIOUS_MODELS_",i,".RData",sep=""))
}
