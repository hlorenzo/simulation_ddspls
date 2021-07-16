source('~/Documents/GitHub/simulation_ddspls/functions.R')
library(ddsPLS2)

# N <- 100
# RESULTSddsPLS <- list()
# ids <-
# for (i in 1:N){
#   cat(i)
#   cat("    ||||    ")
#   data <- get_design_2(seed=i)
#   RESULTSddsPLS[[i]] <- ddsPLS(do.call(cbind,data$Xs),Y = data$Y,n_B = 200,verbose=F,n_lambdas = 40)
#   cat(RESULTSddsPLS[[i]]$R)
#   cat("    ||||    ")
#   cat(round(RESULTSddsPLS[[i]]$Q2,2))
#   cat("\n")
# }

ids <- 1:100
# NCORES <- 10
# idd <- unlist(lapply(1:NCORES,function(ii){rep(ii,10)}))
# `%my_do%` <- `%dopar%`
# cl <- makeCluster(NCORES)
# registerDoParallel(cl)
# RESULTSddsPLS <- foreach(ii=1:NCORES,
#                .packages = "ddsPLS2",
#                .combine='c',
#                .multicombine=TRUE) %my_do% {
#                  models_lambda0 <- list()
#                  iiii <- c(1:100)[which(idd==ii)]
#                  for(i in iiii){
#                    data <- get_design_2(seed=i,sigmaY=0.2)
#                    x <- do.call(cbind,data$Xs)
#                    y <- data$Y
#                    models_lambda0[[i]] <- ddsPLS(x,y,verbose = F,n_B = 200,n_lambdas = 60)
#                  }
#                  models_lambda0
#                }
# stopCluster(cl)

# RESULTSddsPLS2 <- res[which(unlist(lapply(RESULTSddsPLS,length))>0)]
RESULTSddsPLS2 <- list()
for(i in 27:100){
  data <- get_design_2(seed=i,sigmaY=0.2)
  x <- do.call(cbind,data$Xs)
  y <- data$Y
  t0 <- Sys.time()
  RESULTSddsPLS2[[i]] <- ddsPLS(x,y,verbose = F,n_B = 200,n_lambdas = 20)
  po <- do.call(cbind,RESULTSddsPLS2[[i]]$lambda_optim)*1
  po1=po2 <- po; po1[which(po1==0)] <- NA;po2[which(po2==1)] <- NA;po2[which(po2==0)] <- 1;
  bibi <- do.call(cbind,RESULTSddsPLS2[[i]]$results$Q2mean)
  lili <- RESULTSddsPLS2[[i]]$results$lambdas
  matplot(lili,bibi,type="l",lty=3,lwd=3,ylim=range(bibi),main=i)
  matplot(lili,bibi*po1,type="l",lty=1,lwd=3,add=T)
  points(RESULTSddsPLS2[[i]]$lambda,RESULTSddsPLS2[[i]]$Q2,col=1:RESULTSddsPLS2[[i]]$R,pch=16,cex=3)
  cat(i);cat(", R=");cat(RESULTSddsPLS2[[i]]$R);cat(", ");print(Sys.time()-t0); cat("\n")
  print(table(unlist(lapply(RESULTSddsPLS2,function(mo){mo$R}))))
  cat("\n")
}
