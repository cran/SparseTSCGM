
sparse.tscgm <- function (data=data, lam1=lam1, lam2=lam2, B.start = NULL, model.select=c("null","bic"), maxit.out = 100, maxit.in = 500, 
      tol.out = 1e-04, tol.in = 1e-05, tol=1e-05, silent = TRUE) 
{
    model.select = match.arg(model.select)
    nobic = (length(lam1) + length(lam2) == 2)
    doms=(length(lam1)+length(lam2) > 2)

    donull = nobic & (model.select == "null")
    dobic =  doms & (model.select == "bic")
   
    Bhat = NULL
    omega = NULL
    lambda1 = NULL
    lambda2 = NULL
     if (donull) {
        tmp.out = compute.sparse.tscgm(data=data, lam1 = lam1, lam2 = lam2, tol.out = tol.out, 
                 tol.in = tol.in, tol = tol, maxit.out = maxit.out, maxit.in = maxit.in, silent = silent)
        Bhat = tmp.out$Bhat
        omega = tmp.out$omega
    }
   
    else if (dobic) {
        tmp.out = sparse.tscgm.bic(data=data, lam1 = lam1, lam2 = lam2, maxit.out = maxit.out,
                 maxit.in = maxit.in, tol.out = tol.out, tol.in = tol.in, tol=tol,silent = silent)
        Bhat = tmp.out$Bhat
        omega = tmp.out$omega
        lambda1 = tmp.out$lambda1
        lambda2 = tmp.out$lambda2
       
    }
  return(list(Bhat = Bhat, omega = omega, lambda1 = lambda1, 
        lambda2 = lambda2))
}

