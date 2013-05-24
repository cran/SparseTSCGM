\name{sparse.tscgm}
\alias{sparse.tscgm}
\title{
Sparse time series chain graphical models
}
\description{
Computes sparse autoregressive coefficients matrix and precision matrix estimates 
  for time series chain graphical models using SCAD penalty. In time series chain
  graphical models \code{directed} edges or links  are identified by nonzero 
  entries of the autoregressive coefficients matrix and \code{undirected} edges 
  or links  are identified by nonzero entries of the precision matrix.
}
\usage{
sparse.tscgm(data=data, lam1 = lam1, lam2 = lam2, B.start = NULL, 
     model.select = c("null", "bic"), maxit.out = 100,  maxit.in = 500, 
     tol.out = 1e-04, tol.in = 1e-05, tol = 1e-05, silent = TRUE)
}
\arguments{
   \item{data}{Longitudinal data format.
}
  \item{lam1}{
Either a scalar or a vector of tuning parameter values for the 
  SCAD penalty on the precision matrix. If \code{lam1} is a vector with length 
  at least 2, \code{minimum BIC} will be performed to select the optimal tuning
  parameter.
}
  \item{lam2}{
Either a scalar or a vector of tuning parameter values for the 
  SCAD penalty on the autoregression matrix. If \code{lam2} is a vector with length 
  at least 2, \code{minimum BIC} will be performed to select the optimal tuning
  parameter.
}
  \item{B.start}{
The user can optionally specify a starting value for the 
  autoregresive coefficients matrix.
}
  \item{model.select}{
\code{model.select="null"} computes the time series chain graphical model solutions 
  using specified scalar values of the tuning parameters \code{lam1} and \code{lam2},
  \code{model.select="bic"} first selects the optimal SCAD tuning parameter values 
  using the minimum BIC criterion from specified vectors of tuning parameters for 
  \code{lam1} and \code{lam2} then computes the time series chain graphical model   
  solutions based on the optimal tuning parameter values.
  The default is \code{model.select="null".}
}
 \item{maxit.out}{The maximum allowable outer iterations of the algorithm. The
  outer-loop iterates alternatively to estimate the precision matrix and the 
  autoregression matrix until convergence.  
}
  \item{maxit.in}{The maximum allowable iterations of the algorithm that minimizes 
  with respect to the autoregression coefficient matrix.
}
  \item{tol.out}{Convergence tolerance for outer-loop of the algorithm that alternate 
  between the estimation of precision matrix and autoregressive coefficients
  matrix.
}
  \item{tol.in}{Convergence tolerance for inner-loop of the algorithm  that minimizes 
  with respect to the autoregression coefficient matrix.
}
  \item{tol}{Convergence tolerance for the algorithm  that minimizes 
  with respect to the autoregressive coefficients matrix.
}
  \item{silent}{Either \code{TRUE} or \code{FALSE}.  If \code{silent=FALSE} 
  this function will display progress updates to the screen.
}
}
\details{
For description of the objective functions see Abegaz and Wit (2012). 
}
\value{
A list containing:
\item{Bhat}{Autoregressive coefficients matrix estimate. The nonzero estimates of 
the autoregression matrix can be used to represent \code{directed} edges or links  for
chain graphical modeling.}
\item{omega}{Precision matrix estimate. The nonzero estimates of 
the precision matrix can be used to represent \code{undirected} edges or links  for
chain graphical modeling.}
\item{lambda1}{The optimal tuinng parameter for the SCAD penalty on the precision
matrix selected with minimum BIC criterion. If tuning parameter selection 
using BIC is not performed, \code{lambda1} will be \code{NULL}}.
\item{lambda2}{The optimal tuinng parameter for the SCAD penalty on the autoregressive
coefficients matrix selected with minimum BIC criterion. If tuning parameter selection 
using BIC is not performed, \code{lambda2} will be \code{NULL}.} 
}
\references{
Fentaw Abegaz and Ernst wit (2012) Sparse time series chain graphical models  
  for reconstructing genetic networks. Biostatistics, doi:10.1093/biostatistics/kxt005.
  
Rothman, A.J., Levina, E., and Zhu, J. (2010) Sparse multivariate regression 
  with covariance estimation. Journal of Computational and Graphical 
  Statistics. 19: 947--962.
}
\author{
  Fentaw Abegaz and Ernst Wit
}
\note{
For some tuning parameter values the  algorithm may run very slowly. 
 Use \code{silent=FALSE} to check if the algorithm
 is converging before the total iterations exceeds \code{maxit.out}. 
}


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
