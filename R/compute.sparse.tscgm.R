
compute.sparse.tscgm <- function(data=data, lam1=lam1, lam2=lam2, model.select = c("null", "bic"), maxit.out = 100, maxit.in = 500, 
      tol.out = 1e-04, tol.in = 1e-05,  tol = 1e-05,silent=FALSE)
{
  if(is.longitudinal(data)==TRUE)
  {
   time = get.time.repeats(data)$repeats[[1]]
   p_num = dim(data)[2] 
   n_obs = dim(data)[1]/time
   xy.data <-array(NA, c(time,p_num,n_obs))
   data <- as.matrix(data)
   for(i in 1:n_obs){
     for(t in 1:time){
       cc <- 1+(t-1)*n_obs + (i-1)
       xy.data[t,,i] <- data[cc,]
     }
   }
 }
 else cat("Data format is not longitudinal.", "\n")

  X <- xy.data[1:time-1,,]
  Y <- xy.data[2:time,,]
  T <- dim(Y)[1]
  p <- dim(X)[2]
  n <- dim(Y)[3]
  q <- dim(Y)[2] 
  xtyi <- array(NA, c(p,q,n))
  xtxi <- array(NA, c(p,p,n))
  ytyi <- array(NA, c(q,q,n))
      
  for(i in 1:n){
      XX <- X[,,i]
      YY <- Y[,,i]
      xtyi[,,i]=crossprod(XX,YY)
      xtxi[,,i]=crossprod(XX)
      ytyi[,,i]=crossprod(YY)

    }
  xty=apply(xtyi, c(1,2), sum)
  xtx=apply(xtxi, c(1,2), sum)
  yty=apply(ytyi, c(1,2), sum)
     xtxt=apply(xtxi, c(1,2), sum)/(n*T)
  nlam=(n*T)*lam2
  #starting value for matrix autoregressive coefficients	
  old.B = qr.solve(xtx + nlam*diag(p), xty)
  k=0
  mab =  sum(sum(abs(old.B)))
  while(1)
  {
    k=k+1
        ##Initial estimates for the precision matrix using Glasso
        samp.cov = ( yty - t(xty) %*% old.B - t(old.B) %*% xty + t(old.B) %*% xtx %*% old.B)/(n*T)
        g.out=glasso(s=samp.cov, rho=lam1, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)
        old1.om=g.out$wi
	  old1.om.i=g.out$w 
          
        # Calculating SCAD penalty
        wt <- matrix(NA, nrow=q,ncol=q)
        a <- 3.7
        for(i in 1:q){
          for(j in 1:q){
           if(abs(old1.om[i,j]) <= lam1 ) wt[i,j] <- lam1
             else if((lam1 <= abs(old1.om[i,j])) & (abs(old1.om[i,j]) <  a*lam1 )) wt[i,j] <- ((a*lam1-abs(old1.om[i,j]))/((a-1)))
          	   else wt[i,j] <- 0 
        }}
        ##Estimating the precision matrix using Glasso based on SCAD penalty
        g.out1=glasso(s=samp.cov, rho=wt, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE, start="warm", w.init=old1.om.i, wi.init=old1.om)
	    old.om=g.out1$wi
	    old.om.i=g.out1$w 
 
        ##Estimating the auotregressive coeficients based on SCAD penalty
        xtyom=(xty%*%old.om)/(n*T)
        wt1 <- matrix(NA, nrow=p,ncol=q)
        a <- 3.7
        for(i in 1:p){
          for(j in 1:q){
            if(abs(old.B[i,j]) <= lam2 ) wt1[i,j] <- lam2
              else if((lam2 < abs(old.B[i,j])) & (abs(old.B[i,j]) <  a*lam2 )) wt1[i,j] <- ((a*lam2-abs(old.B[i,j]))/((a-1)))
                else wt1[i,j] <- 0 
          }}
        rho2 <- wt1
        warmstart=1
 	  if(k == 1) warmstart=0
 	  B=rblasso(s=xtxt, m=xtyom, om=old.om, nlam=rho2, tol=tol, sbols=mab, maxit=maxit.in, warm=warmstart, B0=old.B)		
  	  bdist = sum(sum(abs(B-old.B)))
   	  old.B=B
          if( (bdist < tol.out*mab) | (k > maxit.out))
    	      break	

   }
  if(silent ==FALSE) cat("Total outer iterations for tscgm : ", k, "\n")
  return(list(Bhat=old.B, omega=old.om))
}

