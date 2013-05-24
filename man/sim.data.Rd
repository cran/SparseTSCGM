\name{sim.data}
\alias{sim.data}
\title{
Simulation from first order autoregressive time series chain graphical models
}
\description{
This function provides repeated multivariate time series simulated data from
random chain graph represented by precision and autoregressive coefficients matrices.
}
\usage{
sim.data(time=time,n_obs=n_obs, n_var=n_var,prob0=prob0)
}
\arguments{
   \item{time}{Number of time points.
}
  \item{n_obs}{Number of repeated observations at each time point.
 }
  \item{n_var}{ Number of variables.
}
  \item{prob0}{The probability that a pair of nodes has an edge. 
 }
}
\value{
A list containing:
\item{data1}{Simulated data in longitudinal data format.}
\item{true_omega}{Precision matrix.}
\item{true_gamma}{Autoregressive coefficients matrix.}
}
\references{
Fentaw Abegaz and Ernst wit (2012) Sparse time series chain graphical models
  for reconstructing genetic networks. Biostatistics, doi:10.1093/biostatistics/kxt005.

Li H., Zhao T., Wang L. Yuan X. and Lui H. (2013), flare: Family of Lasso Regression (R
  package), http://flare.r-forge.r-project.org
}
\author{
  Fentaw Abegaz and Ernst Wit
}

\examples{
#Simulated data from random graph
set.seed(1)
simu <- sim.data(time=10,n_obs=5, n_var=20,prob0=0.15)
Sdata <- simu$data1
prec.mat <- simu$true_omega
autoreg.mat <- simu$true_gamma
}