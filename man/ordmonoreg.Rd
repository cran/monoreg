\name{ordmonoreg}
\alias{ordmonoreg}
\title{Bayesian monotonic regression}
\description{
This function implements the non-parametric Bayesian monotonic regression procedure for 
ordinal outcomes described in Saarela, Rohrbeck & Arjas (2022). 
}
\usage{

ordmonoreg(niter=15000, burnin=5000, adapt=5000, refresh=10, thin=20, 
		   birthdeath=1, logit=FALSE, gam=FALSE, seed=1, rhoa=0.1, rhob=0.1, 
		   deltai=0.5, dlower=0, dupper=1, invprob=1.0, dc=0.0, 
		   predict, include, outcome, axes, covariates=NULL, 
		   cluster=NULL, ncluster=NULL, settozero)
}
\arguments{
\item{niter}{Total number of MCMC iterations.}
\item{burnin}{Number of iterations used for burn-in period.}
\item{adapt}{Number of iterations used for adapting Metropolis-Hastings proposals.}
\item{refresh}{Interval for producing summary output of the state of the MCMC sampler.}
\item{thin}{Interval for saving the state of the MCMC sampler.}
\item{birthdeath}{Number of birth-death proposals attempted at each iteration.}
\item{logit}{Indicator for fitting the model on logit scale.}
\item{gam}{Indicator for fitting non-monotonic generalized additive models for covariate effects.}
\item{seed}{Seed for the random generator.}
\item{rhoa}{Shape parameter of a Gamma hyperprior for the Poisson process rate parameters.}
\item{rhob}{Scale parameter of a Gamma hyperprior for the Poisson process rate parameters.}
\item{deltai}{Range for a uniform proposal for the function level parameters.}
\item{dlower}{Lower bound for the allowed range for the monotonic function.}
\item{dupper}{Upper bound for the allowed range for the monotonic function.}
\item{invprob}{Probability with which to propose keeping the original direction of monotonicity.}
\item{dc}{First parameter of the conditional beta prior to counter spiking behaviour near the origin - 1.}
\item{predict}{Indicator vector for the observations for which absolute risks are calculated (not included in the likelihood expression).}
\item{include}{Indicator vector for  the observations to be included in the likelihood expression.}
\item{outcome}{Vector of outcome variable values.}
\item{axes}{A matrix where the columns specify the covariate axes in the non-parametrically specified regression functions. Each variable here must be scaled to zero-one interval.}
\item{covariates}{A matrix of additional covariates to be included in the linear predictor of the events of interest (optional).}
\item{cluster}{A vector of indicators for cluster membership, numbered from 0 to ncluster-1 (optional).}
\item{ncluster}{Number of clusters (optional).}
\item{settozero}{A zero-one matrix specifying the point process construction. Each row represents a point process,
while columns correspond to the columns of the argument \code{axes}, indicating whether the column is one of the dimensions specifying the domain of the point process.
(See function \code{\link{getcmat}}.)}
}
\value{
A list with elements
\item{steptotal}{A sample of total number of points in the marked point process construction.}
\item{steps}{A sample of the number of points used per each additive component.}
\item{rho}{A sample of the Poisson process rate parameters (one per each point process specified).}
\item{loglik}{A sample of log-likelihood values.}
\item{tau}{A sample of random intercept variances.}
\item{alpha}{A sample of random intercepts.}
\item{beta}{A sample of regression coefficients for the variables specified in the argument \code{covariates}.}
\item{lambda}{A sample of non-parametric regression function levels for each observation specified in the argument \code{predict}.}
\item{pred}{A sample of predicted means for each observation specified in the argument \code{predict}.}
\item{sigmasq}{A sample of autoregressive prior variances.}
}
\references{
Saarela O., Rohrbeck C., Arjas E. (2022). Bayesian non-parametric ordinal regression under a monotonicity constraint. Bayesian Analysis, accepted for publication. arXiv:2007.01390
}
\examples{
library(monoreg)
expit <- function(x) {1/(1+exp(-x))}
logit <- function(p) {log(p)-log(1-p)}
set.seed(1)
nobs <- 500
x <- sort(runif(nobs))
ngrid <- 100
xgrid <- seq(1/ngrid, (ngrid-1)/ngrid, by=1/ngrid)
ngrid <- length(xgrid)
ncat <- 4
beta <- 0.75
disc <- c(Inf, Inf, 0.75, 0.5)
gamma <- c(0,0,0.25,0.5)

surv <- matrix(NA, nobs, ncat)
cdf <- matrix(NA, nobs, ncat)
cols <- c('black','red','blue','green')
for (i in 1:ncat) {
    surv[,i] <- expit(logit((ncat-(i-1))/ncat) + beta * x + 
	                  gamma[i] * (x > disc[i]) - gamma[i] * (x < disc[i]))
    if (i==1)
        plot(x, surv[,i], type='l', col=cols[i], ylim=c(0,1), lwd=2, ylab='S')
    else
        lines(x, surv[,i], col=cols[i], lwd=2)
}
head(surv)

for (i in 1:ncat) {
    if (i<ncat)
        cdf[,i] <- 1.0 - surv[,i+1]
    else
        cdf[,i] <- 1.0
    if (0) {
        if (i==1)
            plot(x, cdf[,i], type='l', col=cols[i], ylim=c(0,1), lwd=2, ylab='F')
        else
            lines(x, cdf[,i], col=cols[i], lwd=2)
    }
}
head(cdf)

u <- runif(nobs)
y <- rep(NA, nobs)
for (i in 1:nobs)
    y[i] <- findInterval(u[i], cdf[i,]) + 1
table(y)

xwindow <- 0.1/2
mw <- matrix(NA, nobs, ncat)
grid <- sort(x)
for (i in 1:nobs) {
    idx <- (x > grid[i] - xwindow) & (x < grid[i] + xwindow)
    for (j in 1:ncat)
        mw[i,j] <- sum(y[idx] >= j)/sum(idx)
}
for (j in 1:ncat) {
    lines(grid, mw[,j], lty='dashed', col=cols[j])
}

# results <- ordmonoreg(niter=15000, burnin=5000, adapt=5000, refresh=10, thin=5, 
results <- ordmonoreg(niter=3000, burnin=1000, adapt=1000, refresh=10, thin=4, 
           birthdeath=1, logit=FALSE, gam=FALSE, seed=1, rhoa=0.1, rhob=0.1, 
           deltai=0.2, dlower=0, dupper=1, invprob=1.0, dc=0.0, 
           predict=c(rep(0,nobs),rep(1,ngrid)), 
           include=c(rep(1,nobs),rep(0,ngrid)), 
           outcome=c(y, rep(1,ngrid)), 
           axes=c(x, xgrid), covariates=NULL, cluster=NULL, ncluster=NULL, 
           settozero=getcmat(1))
		   
s <- results$lambda
dim(s)
lines(xgrid, colMeans(subset(s, s[,1]==0)[,2:(ngrid+1)]))
lines(xgrid, colMeans(subset(s, s[,1]==1)[,2:(ngrid+1)]), col='red')
lines(xgrid, colMeans(subset(s, s[,1]==2)[,2:(ngrid+1)]), col='blue')
lines(xgrid, colMeans(subset(s, s[,1]==3)[,2:(ngrid+1)]), col='green')
legend('bottomright', legend=c(expression(P(Y>=1)), expression(P(Y>=2)), 
expression(P(Y>=3)), expression(P(Y>=4))), 
lwd=2, col=cols)
}
\author{Olli Saarela <olli.saarela@utoronto.ca>, Christian Rohrbeck <cr777@bath.ac.uk>}
