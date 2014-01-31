
sim.data <-
  function(model=c("ar1","ar2"),time=time,n.obs=n.obs, n.var=n.var,prob0=prob0,network=c("random","cluster","hub")){
  model=match.arg(model)
  network=match.arg(network)
  t=time
  n=n.obs
  d=n.var
  if(model=="ar1") {
    if(network=="random") {
      L = flare.tiger.generator(n=n,d=d,graph="random", prob=prob0, seed=1234, vis = FALSE)
      LL = flare.tiger.generator(n=n,d=d,graph="random", prob=prob0, seed=4567, vis = FALSE)
    }
    else if(network=="cluster") {
      L = flare.tiger.generator(n=n,d=d,graph="cluster", prob=prob0, seed=1234, vis = FALSE)
      LL = flare.tiger.generator(n=n,d=d,graph="cluster", prob=prob0, seed=4567, vis = FALSE)
    }
    else if(network=="hub") {
     L = flare.tiger.generator(n=n,d=d,graph="hub", prob=prob0, seed=1234, vis = FALSE)
    LL = flare.tiger.generator(n=n,d=d,graph="hub", prob=prob0, seed=4567, vis = FALSE)
   }
  true_theta = as.matrix(L$omega*L$theta)
  diag(true_theta)=1     ##theta is the precision matrix
  sigma1 <- L$sigma
  mu <- rep(0,d)
  thetaL = as.matrix(LL$omega*LL$theta)
  lwt = 0*lower.tri(thetaL, diag =FALSE)
  upt = thetaL*(1*upper.tri(thetaL, diag = FALSE))
  B=upt+lwt
  ua=rbinom(d,1,0.02)
  uu=runif(d,0,1)
  uau = ua*uu
  diag(B)= uau
  true_gamma=B         #Gamma is the autoregressive coefficient matrix
  ##data generation
  xtn <-array(NA,c(t,d,n))
  xtt <- array(NA,c(t,d,1))
  for(i in 1:n){
    x0 <- rmvnorm(1,mu,sigma1)
    for(j in 1:t){
      et <- rmvnorm(1,mu,sigma1)
      xt <-  x0 %*% B + et
      xtt[j,,] <- xt
      x0 <- xt
   }
   xtn[,,i] <- round(xtt,3)
 }
 xy=matrix(aperm(xtn, c(3,1,2)), ncol=d)
 data1 <- as.longitudinal(xy, repeats=n)
 return(list(data1=data1,theta=true_theta, gamma=true_gamma))
 }
 if(model=="ar2") {
  if(network=="random") {
    L = flare.tiger.generator(n=n,d=d,graph="random", prob=prob0, seed=12346, vis = FALSE)
    LL = flare.tiger.generator(n=n,d=d,graph="random", prob=prob0, seed=45678, vis = FALSE)
    LLL = flare.tiger.generator(n=n,d=d,graph="random", prob=prob0, seed=43219, vis = FALSE)
  }
  else if(network=="cluster") {
   L = flare.tiger.generator(n=n,d=d,graph="cluster", prob=prob0, seed=12346, vis = FALSE)
   LL = flare.tiger.generator(n=n,d=d,graph="cluster", prob=prob0, seed=45678, vis = FALSE)
   LLL = flare.tiger.generator(n=n,d=d,graph="cluster", prob=prob0, seed=43219, vis = FALSE)
  }
 else if(network=="hub") {
  L = flare.tiger.generator(n=n,d=d,graph="hub", prob=prob0, seed=12346, vis = FALSE)
  LL = flare.tiger.generator(n=n,d=d,graph="hub", prob=prob0, seed=45678, vis = FALSE)
  LLL = flare.tiger.generator(n=n,d=d,graph="hub", prob=prob0, seed=43219, vis = FALSE)
 }
 true_theta = as.matrix(L$omega*L$theta)
 diag(true_theta)=1     ##omega is the precision matrix
 sigma1 <- L$sigma
 mu <- rep(0,d)
 omegaL = as.matrix(LL$omega*LL$theta)
 lwt = 0*lower.tri(omegaL, diag =FALSE)
 upt = omegaL*(1*upper.tri(omegaL, diag = FALSE))
 B1=upt+lwt
 ua=rbinom(d,1,0.02)
 uu=runif(d,0,1)
 uau = ua*uu
 diag(B1)= uau
 omegaLL = as.matrix(LLL$omega*LLL$theta)
 lwt0 = 0*lower.tri(omegaLL, diag =FALSE)
 upt0 = omegaLL*(1*upper.tri(omegaLL, diag = FALSE))
 B2=upt0+lwt0
 ua=rbinom(d,1,0.01)
 uu=runif(d,0,1)
 uau = ua*uu
 diag(B2)= uau
#B1and B2 are the autoregressive coefficient matrices
##data generation
 xtn <-array(NA,c(t,d,n))
 xtt <- array(NA,c(t,d,1))
 for(i in 1:n){
  x0 <- rmvnorm(1,mu,sigma1)
  x1 <- rmvnorm(1,mu,sigma1)
  for(j in 1:t){
    et <- rmvnorm(1,mu,sigma1)
    xt <-   x1 %*% B2 +  x0 %*% B1 + et
    xtt[j,,] <- xt
    x0 <- x1
    x1 <- xt
  }
   xtn[,,i] <- round(xtt,3)
  }
  true_gamma <- rbind(B2, B1)
  xy=matrix(aperm(xtn, c(3,1,2)), ncol=d)
  data1 <- as.longitudinal(xy, repeats=n)
  
 return(list(data1=data1,theta=true_theta, gamma=true_gamma))
 }
}

