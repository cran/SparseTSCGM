\name{SparseTSCGM-package}
\alias{SparseTSCGM-package}
\alias{SparseTSCGM}
\docType{package}
\title{
Sparse time series chain graphical models.
}
\description{
Computes sparse autoregressive coefficient and precision matrices for time
 series chain graphical models.
}

\author{
Fentaw Abegaz and Ernst Wit

Maintainer: Fentaw Abegaz <f.abegaz.yazew@rug.nl>
}
\references{
Fentaw Abegaz and Ernst wit (2013) Sparse time series chain graphical models  
  for reconstructing genetic networks. Biostatistics, doi:10.1093/biostatistics/kxt005.
  
Rothman, A.J., Levina, E., and Zhu, J. (2010) Sparse multivariate regression 
  with covariance estimation. Journal of Computational and Graphical 
  Statistics. 19: 947--962.
}
\keyword{package}


\examples{
#Simulated data from random graph
set.seed(1)
simu <- sim.data(time=10,n_obs=5, n_var=20,prob0=0.15)
Sdata <- simu$data1

lam1=0.43
lam2=0.27

res.tscgm <- sparse.tscgm(data=Sdata, lam1=lam1, lam2=lam2, B.start = NULL, 
              model.select="null", maxit.out = 100, maxit.in = 500, 
             tol.out = 1e-04, tol.in = 1e-05, silent = TRUE) 
B <- res.tscgm$Bhat
om <- res.tscgm$omega

### Plot chain graph

  ### plot of undirected graph
  
colnames(om) <-  rownames(om) <- colnames(Sdata)
network_precm <- network(om)
plot.network(network_precm,label = network.vertex.names(network_precm),
  usearrows = FALSE, displayisolates = FALSE, pad = 0.01,label.pad = 0.3, 
  label.col = 6, vertex.col = 5,vertex.cex = 1.5,
  edge.col = 4, mode = "fruchtermanreingold", interactive=FALSE)

   ### Plot of directed graph
 
colnames(B) <-  rownames(B) <- colnames(Sdata)
network.autoregm <- network(B)
plot.network(network.autoregm,label = network.vertex.names(network.autoregm),
  usearrows = TRUE, displayisolates = FALSE, pad = 0.01,label.pad = 0.3, 
  label.col = 6, vertex.col = 5,vertex.cex = 1.5,
  edge.col = 3, mode = "fruchtermanreingold", interactive=FALSE)
}
