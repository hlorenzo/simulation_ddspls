get_design_1 <- function(n=50,sqrt_1_minus_sig2=0.99,p=1000,q=3){
  # Structure
  alpha3 <- 1/sqrt(3)
  alpha2 <- 1/sqrt(2)
  repX <- 50
  A1 <- c(rep(alpha3,repX),rep(0,p-repX))
  A2 <- c(rep(0,repX),rep(alpha2,repX),rep(0,p-2*repX))
  A <- matrix(c(rep(A1,3),rep(A2,2)),nrow = 5,byrow = T)*sqrt_1_minus_sig2
  D1 <- c(rep(alpha3,1),rep(0,q-1))
  D2 <- c(rep(0,1),rep(alpha2,1),rep(0,q-2))
  D <- matrix(c(rep(D1,3),rep(D2,2)),nrow = 5,byrow = T)*sqrt_1_minus_sig2
  # Observations
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
