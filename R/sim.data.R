sim.data <-
  function(time=time,n_obs=n_obs, n_var=n_var,prob0=prob0){
 t=time
n=n_obs
d=n_var

L = flare.tiger.generator(n=n,d=d,graph="random", prob=prob0, seed=1234, vis = FALSE)
true_omega = as.matrix(L$omega*L$theta)
diag(true_omega)=1     ##omega is the precision matrix
sigma1 <- L$sigma
mu <- rep(0,d)

LL = flare.tiger.generator(n=n,d=d,graph="random", prob=prob0, seed=4567, vis = FALSE)
plot(LL)
set.seed(2)
omegaL = as.matrix(LL$omega*LL$theta)
lwt = 0*lower.tri(omegaL, diag =FALSE)
upt = omegaL*(1*upper.tri(omegaL, diag = FALSE))
B=upt+lwt
ua=rbinom(d,1,0.2)
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
 return(list(data1=data1,true_omega=true_omega, true_gamma=true_gamma))
}   