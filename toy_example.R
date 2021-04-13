eps1=eps2=eps3=epsY=eps <- 0.95#0.8#

ff <- function(cc){out <- cc/sqrt(sum(cc^2));if(is.na(out[1])) out <- cc;out}
A <- apply(cbind(c(1,1),c(0,1),c(0,0)),2,ff)
d <- nrow(A)
A <- eps1*cbind(A,matrix(0,nrow = d,ncol = 2))
p1  <- 50; p <- 1000
R1 <- 1; R2 <- 0;d <- R1+R2
A <- matrix(rep(c(rep(1,p1),rep(0,p-p1)),R1),nrow = R1,byrow = T)

A <- eps1*A
D <- matrix(rep(c(1),R1),nrow = R1,byrow = T)
D <- eps1*D

library(ddsPLS2)

ns <- c(50,100,200)
LAMBDAS <- seq(0,1,length.out = 66)
ls <- matrix(LAMBDAS,ncol = 1)
res_s = Xs_s = Y_s <- list()

L_total <- ncol(D)+ncol(A)
ii<-1
for(ii in 1:3){
  n <- ns[ii]
  psi <- MASS::mvrnorm(n,mu = rep(0,d+L_total),Sigma = diag(d+L_total))
  phi <- psi[,1:d,drop=F];pt <- d
  SIs <- lapply(list(A,D),function(M){do.call(cbind,lapply(sqrt(1-diag(crossprod(M))),function(sisi){rep(sisi,n)}))})
  Xs <- phi%*%A + SIs[[1]]*psi[,pt+1:ncol(A),drop=F];pt <- pt + ncol(A)
  Y <- phi%*%D + SIs[[2]]*psi[,pt+1:ncol(D),drop=F];pt <- pt + ncol(Y)
  Y_s[[ii]] <- Y
  Xs_s[[ii]] <- Xs
  res_s[[ii]] <- ddsPLS(Xs,Y,lambdas=LAMBDAS,n_B=50,NCORES=1,verbose=T)
}

nb_variables_s <- list()
plot(-10,-10,xlim = c(0,66),ylim = c(0,1000))
for(ii in 2:3){
  nb_variables_s[[ii]] <- rep(0,66)
  Xs_s[[ii]] -> Xs
  Y_s[[ii]] -> Y
  for(i in 1:length(LAMBDAS)){
    lam <- LAMBDAS[i]
    res_i <- ddsPLS(Xs,Y,lambdas=lam,n_B=50,NCORES=1,verbose=F)
    if(!is.null(res_i$model)){
      nb_variables_s[[ii]][i] <- sum((rowSums(abs(res_i$model$U))>1e-9)*1)
    }else{
      nb_variables_s[[ii]][i] <- 0
    }
    points(i,nb_variables_s[[ii]][i])
  }
}

######################################################################
postscript("/Users/hlorenzo/Documents/GitHub/simu_facile.eps",
           width=10, height=3.5, onefile=TRUE, horizontal=F)
par(layout(matrix(c(1,2,3,rep(4,2),5),nrow = 2,byrow = T)))
par(mar=c(4,2.9,2,2),cex=0.5);cols <- c(1,1,1);cols_all <- c(4,6,2)
for(ii in 1:3){
  res_i <- res_s[[ii]]
  id_i <- 1:length(LAMBDAS)
  id_i_plo <- seq(5,55,by = 10)-2
  id_max_i_q2 <- which.max(res_i$results$Q2mean[[1]])
  max_i_q2 <- max(res_i$results$Q2mean[[1]])
  col_i <- cols[ii]
  matplot(LAMBDAS[id_i],cbind(res_i$results$R2mean[[1]],
                              res_i$results$Q2mean[[1]])[id_i,],
          main=bquote(n==.(ns[ii])),lwd=1,xlab=expression(lambda),
          ylab="",type="l",lty=c(1,2),ylim=c(0.8,0.96),col=col_i)
  points(LAMBDAS[id_max_i_q2],max_i_q2,pch=11,cex=2,col = cols_all[ii])
  points(c(res_i$lambda[1],res_i$lambda[1]),
         c(res_i$results$R2mean[[1]][which(LAMBDAS==res_i$lambda[[1]])],
           res_i$results$Q2mean[[1]][which(LAMBDAS==res_i$lambda[[1]])]),type="l",col = cols_all[ii])
  points(c(res_i$lambda[1],res_i$lambda[1]),
         c(0,
           res_i$results$Q2hmean[[1]][which(LAMBDAS==res_i$lambda[1])]),type="l",lty=3,col = cols_all[ii])
  text(x=res_i$lambda[1],y = 0.85,col = cols_all[ii],
       labels = bquote(hat(lambda)[1]==.(round(res_i$lambda[1],2) )))
  abline(h=0.9025,lwd=0.8,col="gray30")
  text(0.6,0.91 ,expression(gamma^"*"~"="~1-epsilon^2~"="~0.9025),col="gray30")
  legend("top",ncol = 3,bty="n",
         lty=c(1,2,NA),
         lwd=c(1,1,NA),
         legend = c(expression(bar("R")["B,1"]^"2"),
                    expression(bar("Q")["B,1"]^"2"),
                    expression("Max of"~Q^2)),
         pch=c(NA,NA,11),cex=c(rep(0.9,2),0.9))
}


######################################################################
res_50<- res_s[[1]];nb_variables_50 <- nb_variables_s[[1]];
id_max_50_q2 <- which.max(res_50$results$Q2mean[[1]]);
max_50_q2 <- max(res_50$results$Q2mean[[1]])
id_50_plo <- c(seq(5,64,by = 5)-2,60)
res_100<- res_s[[2]];nb_variables_100 <- nb_variables_s[[2]];  id_max_100_q2 <- which.max(res_100$results$Q2mean[[1]]);max_100_q2 <- max(res_100$results$Q2mean[[1]])
id_100_plo <- c(seq(5,64,by = 5),60)
res_200<- res_s[[3]];nb_variables_200 <- nb_variables_s[[3]];  id_max_200_q2 <- which.max(res_200$results$Q2mean[[1]]);max_200_q2 <- max(res_200$results$Q2mean[[1]])
id_200_plo <- c(seq(5,64,by = 5)+2,62)
plot(-10,-10,xlim=c(0,1),ylim=c(0,sqrt(1000)),
     xlab=expression(lambda),
     yaxt="n",
     ylab="",main="Number of selected variables (quadratic scale)")
abline(h=sqrt(50),lty=1)
abline(h=0,lty=2)
axis(2,at = sqrt(seq(0,1000,by = 200)),gap.axis = 0.01,labels = seq(0,1000,by = 200),las=2)
axis(2,at = sqrt(50),labels = expression(bold(50)),las=2)
matplot(cbind(ls[id_200_plo,1],ls[id_100_plo,1],ls[id_50_plo,1]),
        sqrt(cbind(nb_variables_200[id_200_plo],nb_variables_100[id_100_plo],nb_variables_50[id_50_plo])),
        type="p",pch=c(4,2,1),col=c(2,6,4),add=T,cex=0.7)
matplot(ls[,1],sqrt(cbind(nb_variables_200,nb_variables_100,nb_variables_50)),
        type="l",col=c(2,6,4),add=T,lty=1)
points(LAMBDAS[id_max_50_q2],0,pch=11,cex=2,col=4)
points(LAMBDAS[id_max_100_q2],0,pch=11,cex=2,col=6)
points(LAMBDAS[id_max_200_q2],0,pch=11,cex=2,col=2)
abline(v=c(res_50$optimal_parameters$lambdas,
           res_100$optimal_parameters$lambdas,
           res_200$optimal_parameters$lambdas),col=c(4,6,2),lty=2,lwd=1/2)
text(x=res_200$lambda,y = 25,col=2,
     labels = bquote(hat(lambda)[1](n==200)==.(round(res_200$lambda,2) )))
text(x=res_100$lambda,y = 20,col=6,
     labels = bquote(hat(lambda)[1](n==100)==.(round(res_100$lambda,2) )))
text(x=res_50$lambda,y = 15,col=4,
     labels = bquote(hat(lambda)[1](n==50)==.(round(res_50$lambda,2) )))
legend("topright",ncol = 3,bty="n",
       col=c(2,6,4),
       border=c(1,2,4),
       legend = c("n=200","n=100","n=50"),
       pch=c(4,2,1))

########################################################################
b50<-res_50$model$B[1:50,1]
b100<-res_100$model$B[1:50,1]
b200<-res_200$model$B[1:50,1]
matou <- cbind(b50,b100,b200)
colnames(matou) <- c("n=50","n=100","n=200")
boxplot(matou,main="50 first regression coefficients",col=c(4,6,2),ylab="",yaxt="n")
axis(2,(19:23)/1000,(19:23)/1000,las=2)
abline(h=1/50)

dev.off()






######################################################################
postscript("/Users/hlorenzo/Documents/GitHub/simu_facile_fr.eps",
           width=10, height=3.5, onefile=TRUE, horizontal=F)
par(layout(matrix(c(1,2,3,rep(4,2),5),nrow = 2,byrow = T)))
par(mar=c(4,2.9,2,2),cex=0.5);cols <- c(1,1,1);cols_all <- c(4,6,2)
for(ii in 1:3){
  res_i <- res_s[[ii]]
  id_i <- res_i$id_ALL_TEST_h[[1]]
  id_i_plo <- seq(5,55,by = 10)-2
  id_max_i_q2 <- which.max(res_i$bootstrap$Q2_h_star[[1]])
  max_i_q2 <- max(res_i$bootstrap$Q2_h_star[[1]])
  col_i <- cols[ii]
  matplot(LAMBDAS[id_i],cbind(res_i$bootstrap$vars_h_boot[[1]],
                              res_i$bootstrap$Q2_h_star[[1]])[id_i,],
          main=bquote(n==.(ns[ii])),lwd=1,xlab=expression(lambda),
          ylab="",type="l",lty=c(1,2),ylim=c(0.8,0.96),col=col_i)
  points(LAMBDAS[id_max_i_q2],max_i_q2,pch=11,cex=2,col = cols_all[ii])
  points(c(res_i$optimal_parameters$lambdas,res_i$optimal_parameters$lambdas),
         c(res_i$bootstrap$Q2_h_star[[1]][which(LAMBDAS==res_i$optimal_parameters$lambdas)],
           res_i$bootstrap$vars_h_boot[[1]][which(LAMBDAS==res_i$optimal_parameters$lambdas)]),type="l",col = cols_all[ii])
  points(c(res_i$optimal_parameters$lambdas,res_i$optimal_parameters$lambdas),
         c(0,
           res_i$bootstrap$Q2_h_star[[1]][which(LAMBDAS==res_i$optimal_parameters$lambdas)]),type="l",lty=3,col = cols_all[ii])
  text(x=res_i$optimal_parameters$lambdas,y = 0.85,col = cols_all[ii],
       labels = bquote(hat(lambda)[1]==.(round(res_i$optimal_parameters$lambdas,2) )))
  # abline(h=0.9025,lwd=0.8,col="gray30")
  # text(0.6,0.91 ,expression(gamma^"*"~"="~1-epsilon^2~"="~0.9025),col="gray30")
  legend("top",ncol = 3,bty="n",
         lty=c(1,2,NA),
         lwd=c(1,1,NA),
         legend = c(expression(bar("R")["B,1"]^"2"),
                    expression(bar("Q")["B,1"]^"2"),
                    expression("Max de"~Q^2)),
         pch=c(NA,NA,11),cex=c(rep(0.9,2),0.9))
}


######################################################################
res_50<- res_s[[1]];nb_variables_50 <- nb_variables_s[[1]];  id_max_50_q2 <- which.max(res_50$bootstrap$Q2_h_star[[1]]);max_50_q2 <- max(res_50$bootstrap$Q2_h_star[[1]])
id_50_plo <- c(seq(5,64,by = 5)-2,60)
res_100<- res_s[[2]];nb_variables_100 <- nb_variables_s[[2]];  id_max_100_q2 <- which.max(res_100$bootstrap$Q2_h_star[[1]]);max_100_q2 <- max(res_100$bootstrap$Q2_h_star[[1]])
id_100_plo <- c(seq(5,64,by = 5),60)
res_200<- res_s[[3]];nb_variables_200 <- nb_variables_s[[3]];  id_max_200_q2 <- which.max(res_200$bootstrap$Q2_h_star[[1]]);max_200_q2 <- max(res_200$bootstrap$Q2_h_star[[1]])
id_200_plo <- c(seq(5,64,by = 5)+2,62)
plot(-10,-10,xlim=c(0,1),ylim=c(0,sqrt(1000)),
     xlab=expression(lambda),
     yaxt="n",
     ylab="",main="Nombre de variables sélectionnées (échelle quadratique)")
abline(h=sqrt(50),lty=1)
abline(h=0,lty=2)
axis(2,at = sqrt(seq(0,1000,by = 200)),gap.axis = 0.01,labels = seq(0,1000,by = 200),las=2)
axis(2,at = sqrt(50),labels = expression(bold(50)),las=2)
matplot(cbind(ls[id_200_plo,1],ls[id_100_plo,1],ls[id_50_plo,1]),
        sqrt(cbind(nb_variables_200[id_200_plo],nb_variables_100[id_100_plo],nb_variables_50[id_50_plo])),
        type="p",pch=c(4,2,1),col=c(2,6,4),add=T,cex=0.7)
matplot(ls[,1],sqrt(cbind(nb_variables_200,nb_variables_100,nb_variables_50)),
        type="l",col=c(2,6,4),add=T,lty=1)
points(LAMBDAS[id_max_50_q2],0,pch=11,cex=2,col=4)
points(LAMBDAS[id_max_100_q2],0,pch=11,cex=2,col=6)
points(LAMBDAS[id_max_200_q2],0,pch=11,cex=2,col=2)
abline(v=c(res_50$optimal_parameters$lambdas,
           res_100$optimal_parameters$lambdas,
           res_200$optimal_parameters$lambdas),col=c(4,6,2),lty=2,lwd=1/2)
text(x=res_200$optimal_parameters$lambdas,y = 25,col=2,
     labels = bquote(hat(lambda)[1](n==200)==.(round(res_200$optimal_parameters$lambdas,2) )))
text(x=res_100$optimal_parameters$lambdas,y = 20,col=6,
     labels = bquote(hat(lambda)[1](n==100)==.(round(res_100$optimal_parameters$lambdas,2) )))
text(x=res_50$optimal_parameters$lambdas,y = 15,col=4,
     labels = bquote(hat(lambda)[1](n==50)==.(round(res_50$optimal_parameters$lambdas,2) )))
legend("topright",ncol = 3,bty="n",
       col=c(2,6,4),
       border=c(1,2,4),
       legend = c("n=200","n=100","n=50"),
       pch=c(4,2,1))

########################################################################
b50<-res_50$B_cbind[1:50,1]
b100<-res_100$B_cbind[1:50,1]
b200<-res_200$B_cbind[1:50,1]
matou <- cbind(b50,b100,b200)
colnames(matou) <- c("n=50","n=100","n=200")
boxplot(matou,main="50 premiers coefficients de régression",col=c(4,6,2),ylab="",yaxt="n")
axis(2,(19:23)/1000,(19:23)/1000,las=2)
abline(h=1/50)

dev.off()
