
sparse.tscgm.bic <- 
function(data=data, lam1 = lam1, lam2 = lam2, 
            maxit.out = maxit.out, maxit.in = maxit.in, tol.out = tol.out, 
            tol.in = tol.in, tol=tol, silent = silent)
{
	###Determining the tuning or penalty parameter
	 lam.vec.1 = lam1
	 lam.vec.2 = lam2
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

  X <- xy.data[1:time-1,,]
  Y <- xy.data[2:time,,]     
	lamR =length( lam.vec.1)*length(lam.vec.2)
	BICh = matrix(NA,lamR,1)
	lam1h  = matrix(NA,lamR,1)
	lam2h = matrix(NA,lamR,1)
     
  

	uv <- 0
	for(u in 1:length(lam.vec.1)){
 	 for(v in 1:length(lam.vec.2)){
    	 	uv <- uv+1
     		outscad_s <- compute.sparse.tscgm(data=data, lam1=lam.vec.1[u], lam2=lam.vec.2[v], 
                    maxit.out = maxit.out, maxit.in = maxit.in, tol.out = tol.out, tol.in = tol.in,tol=tol, silent = TRUE) 
     		PO = outscad_s$omega 
     		PB = outscad_s$Bhat 
     		T <- dim(Y)[1]
    		p<- dim(X)[2]
    		n <- dim(Y)[3]
    		 	q<- dim(Y)[2]
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
    
   	 ##log likelihood
    	 lik1 <- log(det(PO))
   	 WS = (yty - t(xty) %*% PB - t(PB) %*% xty + t(PB) %*% xtx %*% PB)/(n*T)
    	 lik2 <- sum(diag(WS%*%PO))
    	 LLk <-  (n*T/2)*(lik1-lik2) 
   
    	###Number of nonzero entries
    	diag(PO)=0
    	pdO = sum(sum(PO !=0)) 
    	pdB = sum(sum(PB !=0)) 

    	###BIC caclulation 
    	BICh[uv,1] <- -2*LLk + log(n*T)*(pdO/2 + q + pdB) #+(pdO/2 + q + pdB)*4*0.5*log(q+p)

    	lam1h[uv,1] <- lam.vec.1[u]
   	lam2h[uv,1] <- lam.vec.2[v]
   }}
   res.scad <- cbind(lam1h, lam2h, BICh)
   ### Determining the minimum Bic and optimal penalty parameter
   bicm <- min(res.scad[,3])
   for(i in 1:lamR){
     if(res.scad[i,3]==bicm){
       lam1.opt <- res.scad[i,1]
       lam2.opt <- res.scad[i,2]
   }}
 
   tmp.out= compute.sparse.tscgm(data=data, lam1=lam1.opt, lam2=lam2.opt,  model.select="null",
	           tol.out=tol.out, tol.in=tol.in, maxit.out=maxit.out, maxit.in=maxit.in, silent=silent) 
   best.B = tmp.out$Bhat
   best.omega = tmp.out$omega
	
   return(list(Bhat=best.B, omega=best.omega, lambda1=lam1.opt, lambda2=lam2.opt))
}

