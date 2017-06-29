\name{monosurv}
\alias{monosurv}
\title{Bayesian monotonic regression for time-to-event outcomes}
\description{
This function implements an extended version of the Bayesian monotonic regression
procedure described in Saarela & Arjas (2011), allowing for multiple additive monotonic
components, and time-to-event outcomes through case-base sampling. Logistic/multinomial
regression is fitted if no time variable is present. The extension and its applications, 
including estimation of absolute risks, are described in Saarela & Arjas (2015).
The example below does logistic regression; for an example of modeling a time-to-event outcome, 
please see the documentation for the dataset \code{\link{risks}}.
}
\usage{
monosurv(niter=15000, burnin=5000, adapt=5000, refresh=10, thin=5, 
         birthdeath=10, timevar=0, seed=1, rhoa=0.1, rhob=0.1, 
         years=NULL, deltai=0.1, drange=2.0, predict, include, 
         casestatus, sprob=NULL, offset=NULL, tstart=NULL, axes, 
         covariates, ccovariates=NULL, settozero, package, cr=NULL)
}
\arguments{
\item{niter}{Total number of MCMC iterations.}
\item{burnin}{Number of iterations used for burn-in period.}
\item{adapt}{Number of iterations used for adapting Metropolis-Hastings proposals.}
\item{refresh}{Interval for producing summary output of the state of the MCMC sampler.}
\item{thin}{Interval for saving the state of the MCMC sampler.}
\item{birthdeath}{Number of birth-death proposals attempted at each iteration.}
\item{timevar}{Number identifying the column in argument \code{axes} representing a time variable. 
Zero if no time variable is present, in which case a logistic/multinomial regression is fitted (instead of a hazard regression).}
\item{seed}{Seed for the random generator.}
\item{rhoa}{Shape parameter of a Gamma hyperprior for the Poisson process rate parameters.}
\item{rhob}{Scale parameter of a Gamma hyperprior for the Poisson process rate parameters.}
\item{years}{Time period over which absolute risks are calculated (on a time scale scaled to zero-one interval; optional).}
\item{deltai}{Range for a uniform proposal for the function level parameters.}
\item{drange}{Allowed range for the monotonic function components, on log-rate/log-odds scale.}
\item{predict}{Indicator vector for the observations for which absolute risks are calculated (not included in the likelihood expression).}
\item{include}{Indicator vector for  the observations to be included in the likelihood expression.}
\item{casestatus}{An integer vector indicating the case status (0=censoring, 1=event of interest, 2=competing event).}
\item{sprob}{Vector of sampling probabilities/rates for each person/person-moment (optional).}
\item{offset}{Vector of offset terms to be included in the linear predictor of the events of interest (optional).}
\item{tstart}{Vector of entry times on a time scale scaled to zero-one interval (optional).}
\item{axes}{A matrix where the columns specify the covariate axes in the non-parametrically specified regression functions. 
Each variable here must be scaled to zero-one interval.}
\item{covariates}{A matrix of additional covariates to be included in the linear predictor of the events of interest. 
Must include at least a vector of ones to specify an intercept term.}
\item{ccovariates}{A matrix of additional covariates to be included in the linear predictor of the competing events (optional).}
\item{settozero}{A zero-one matrix specifying the point process construction. Each row represents a point process,
while columns correspond to the columns of the argument \code{axes}, indicating whether the column is one of the dimensions specifying the domain of the point process.
(See function \code{\link{getcmat}}.)}
\item{package}{An integer vector specifying the additive component into which each point process (row) specified in argument \code{settozero} is placed.}
\item{cr}{A zero-one vector indicating the additive components to be placed in the linear predictor of the competing causes (optional).}
}
\value{
A list with elements
\item{steptotal}{A sample of total number of points in the marked point process construction.}
\item{steps}{A sample of the number of points used per each additive component.}
\item{rho}{A sample of the Poisson process rate parameters (one per each point process specified).}
\item{loglik}{A sample of log-likelihood values.}
\item{beta}{A sample of regression coefficients for the variables specified in the argument \code{covariates}.}
\item{betac}{A sample of regression coefficients for the variables specified in the argument \code{ccovariates}.}
\item{phi}{A sample of non-parametric regression function levels for each observation specified in the argument \code{predict}.}
\item{risk}{A sample of absolute risks of the event of interest for each observation specified in the argument \code{predict}.}
\item{crisk}{A sample of absolute risks of the competing event for each observation specified in the argument \code{predict}.}
}
\references{
Saarela O., Arjas E. (2011). A method for Bayesian monotonic multiple regression. Scandinavian Journal of Statistics, 38:499--513.

Saarela O., Arjas E. (2015). Non-parametric Bayesian hazard regression for chronic disease risk assessment. Scandinavian Journal of Statistics, 42:609--626.
}
\examples{
library(monoreg)
set.seed(1)
# nobs <- 1000
nobs <- 50
x1 <- runif(nobs)
x2 <- runif(nobs)

# 6 different monotonic regression surfaces:
# mu <- sqrt(x1)
mu <- 0.5 * x1 + 0.5 * x2
# mu <- pmin(x1, x2)
# mu <- 0.25 * x1 + 0.25 * x2 + 0.5 * (x1 + x2 > 1.0)
# mu <- 0.25 * x1 + 0.25 * x2 + 0.5 * (pmax(x1, x2) > 0.5)
# mu <- ifelse((x1 - 1.0)^2 + (x2 - 1.0)^2 < 1.0, sqrt(1.0 - (x1 - 1.0)^2 - (x2 - 1.0)^2), 0.0)

y <- rbinom(nobs, 1, mu)

# results <- monosurv(niter=15000, burnin=5000, adapt=5000, refresh=10, 
results <- monosurv(niter=5000, burnin=2500, adapt=2500, refresh=10, 
                    thin=5, birthdeath=10, seed=1, 
                    rhoa=0.1, rhob=0.1, deltai=0.5, drange=10.0, 
                    predict=rep(1.0, nobs), include=rep(1.0, nobs), 
                    casestatus=y, axes=cbind(x1,x2), covariates=rep(1.0, nobs),
                    settozero=getcmat(2), package=rep(1,3))

# pdf(file.path(getwd(), 'pred3d.pdf'), width=6.0, height=6.0, paper='special')
op <- par(mar=c(2,2,0,0), oma=c(0,0,0,0), mgp=c(2.5,1,0), cex=0.75)
pred <- colMeans(results$risk)
idx <- order(pred, decreasing=TRUE)

tr <- persp(z=matrix(c(NA,NA,NA,NA), 2, 2), zlim=c(0,1), 
            xlim=c(0,1), ylim=c(0,1),
            ticktype='detailed', theta=-45, phi=25, ltheta=25, 
            xlab='X1', ylab='X2', zlab='mu')
for (i in 1:nobs) {
    lines(c(trans3d(x1[idx[i]], x2[idx[i]], 0.0, tr)$x, 
            trans3d(x1[idx[i]], x2[idx[i]], pred[idx[i]], tr)$x),
          c(trans3d(x1[idx[i]], x2[idx[i]], 0.0, tr)$y, 
            trans3d(x1[idx[i]], x2[idx[i]], pred[idx[i]], tr)$y),
          col='gray70')
}
points(trans3d(x1[idx], x2[idx], pred[idx], tr), pch=21, bg='white')
par(op)
# dev.off()
}
\author{Olli Saarela <olli.saarela@utoronto.ca>}
