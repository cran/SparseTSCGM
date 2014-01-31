\name{sim.data}
\alias{sim.data}
\title{
Multivariate time series simulation with chain graphical models
}
\description{
Generates sparse vector autoregressive coefficients matrices and
  precision matrix from various network structures and using these matrices
  generates  repeated multivariate  time series dataset.
  }
\usage{
 sim.data(model=c("ar1","ar2"),time=time,n.obs=n.obs, n.var=n.var,
     prob0=prob0, network=c("random","cluster","hub"))
}
\arguments{
  \item{model}{
Specifies the order of vector autoregressive models. Vector autoregressive
  model of order 1 is applied if \code{model} = "ar1" and Vector autoregressive
  model of order 2 is applied if  method = "ar2".
}
  \item{time}{ Number of time points.}
  \item{n.obs}{ Number of observations or replicates.}
  \item{n.var}{Number of variables.}
  \item{prob0}{Initial sparsity level.}
  \item{network}{
Specifies the type of network structure. This could be random, scale-free or hub
structures. Details on simultions from the various network structures can be
found in the R package \code{flare}.
}
  }
\value{
A list containing:
\item{theta}{Sparse precision matrix. }
\item{gamma}{Sparse autoregressive coefficients matrix.}
\item{data1}{Repeated multivariate time series data in longitudinal format.}
}
\author{
  Fentaw Abegaz and Ernst Wit
}

\examples{
set.seed(321)
datas <- sim.data(model="ar1", time=20,n.obs=20, n.var=20,prob0=0.15,
         network="hub")
data.ts <-  datas$data1
prec_true <- datas$theta
autoR_true <- datas$gamma
}
