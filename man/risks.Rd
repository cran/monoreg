\name{risks}
\alias{risks}
\docType{data}
\title{Absolute risks from 7 survival models}
\description{This example dataset includes time-to-event outcomes and
absolute risks for reproducing the ROC curves in Figure 3 of Saarela & Arjas (2015), please see the example code below.}
\usage{data(risks)}
\format{
A data frame containing the elements
\describe{
\item{tstop}{Age at the end of the follow-up (scaled to zero-one interval)}
\item{censvar}{Case status (0=censoring, 1=CVD event, 2=other death).}
\item{tstart}{Age at the start of the follow-up (scaled to zero-one interval).}
\item{model1}{Absolute risk from model 1.}
\item{model2}{Absolute risk from model 2.}
\item{model3}{Absolute risk from model 3.}
\item{model4}{Absolute risk from model 4.}
\item{model5}{Absolute risk from model 5.}
\item{model6}{Absolute risk from model 6.}
\item{model7}{Absolute risk from model 7.}
}
}
\references{
Saarela O., Arjas E. (2015). Non-parametric Bayesian hazard regression for chronic disease risk assessment. Scandinavian Journal of Statistics, 42:609--626.
}
\examples{
\dontrun{
rm(list=ls())
library(monoreg)
library(eha)

# Read the example data:

data(risks)
ftime <- (risks$tstop - risks$tstart)/max(risks$tstop - risks$tstart)
ftime <- ifelse(ftime == 0, ftime + 0.00001, ftime)
censvar <- risks$censvar
case <- censvar == 1
dcase <- censvar == 2
a <- risks$tstart
a2 <- a^2
nobs <- nrow(risks)

# Fit a simple parametric model to remove the effect of baseline age:

wmodel <- weibreg(Surv(ftime, case) ~ a + a2, shape=1)
summary(wmodel)
wcoef <- c(coef(wmodel), 0.0)
wlp <- crossprod(t(cbind(a, a2)), wcoef[1:(length(wcoef)-2)])
wpar <- exp(wcoef[(length(wcoef)-1):length(wcoef)])

whaz <- function(x, bz) {
    return(exp(bz) * (wpar[2] / wpar[1]) * (x / wpar[1])^(wpar[2] - 1))
}
chint <- function(x, bz) {
    return(exp(bz) * (x / wpar[1])^wpar[2])
}
croot <- function(x, bz, c) {
    return(chint(x, bz) - c)
}

# Age-matched case-base sample for model validation:

set.seed(1)
crate <- chint(ftime, wlp)
csrate <- cumsum(crate)
m <- sum(case) * 10
persons <- rep(1:nobs, rmultinom(1, m, crate/sum(crate)))
moments <- rep(NA, m)
for (i in 1:m) {
    u <- runif(1, 0.0, crate[persons[i]])
    moments[i] <- uniroot(croot, c(0.0, ftime[persons[i]]), c=u, 
                          bz=wlp[persons[i]])$root
}
plot(ecdf(risks$tstart[case]), pch=20, col='red')
plot(ecdf(risks$tstart[persons]), pch=20, col='blue', add=TRUE)

rate <- whaz(moments, wlp[persons])
mrate <- mean(rate)

d <- c(rep(0, m), rep(1, sum(censvar == 1)), rep(2, sum(censvar == 2)), censvar)
mom <- c(moments, ftime[censvar == 1], ftime[censvar == 2], rep(1.0, nobs))
per <- c(persons, (1:nobs)[censvar == 1], (1:nobs)[censvar == 2], 1:nobs)

include <- rep(c(1,0), c(m + sum(censvar == 1) + sum(censvar == 2), nobs))
predict <- as.numeric(!include)

offset <- log(sum(crate)/(m * whaz(mom, wlp[per])))
moffset <- rep(log(sum(crate)/(m * mrate)), length(mom))

sprob <- 1/exp(offset)
msprob <- 1/exp(moffset)

stz <- getcmat(2)
settozero <- rbind(stz[1,], stz[1,], stz[2:3,], stz[2:3,])
package <- 1:nrow(settozero)
cr <- c(1,0,rep(1,2),rep(0,2))

# Fit models removing the age effect:

agecir <- matrix(NA, nobs, 7)
for (i in 1:7) {
    agecir[,i] <- as.numeric(colMeans(
monosurv(niter=15000, burnin=5000, adapt=5000, refresh=10, thin=5, 
         birthdeath=10, timevar=1, seed=1, rhoa=0.1, rhob=0.1, 
         years=1.0, deltai=0.1, drange=6.0, predict=predict, include=include, 
         casestatus=d, sprob=msprob, offset=NULL, tstart=NULL, 
         axes=cbind(mom, risks[per,paste('model', i, sep='')]), 
         covariates=rep(1.0, length(per)), ccovariates=rep(1.0, length(per)), 
         settozero=settozero, package=package, cr=cr)$risk))
    print(i)
}

# Fit models without removing the age effect:

cir <- matrix(NA, nobs, 7)
for (i in 1:7) {
    cir[,i] <- as.numeric(colMeans(
monosurv(niter=15000, burnin=5000, adapt=5000, refresh=10, thin=5, 
         birthdeath=10, timevar=1, seed=1, rhoa=0.1, rhob=0.1, 
         years=1.0, deltai=0.1, drange=6.0, predict=predict, include=include, 
         casestatus=d, sprob=sprob, offset=NULL, tstart=NULL, 
         axes=cbind(mom, risks[per,paste('model', i, sep='')]), 
         covariates=rep(1.0, length(per)), ccovariates=rep(1.0, length(per)), 
         settozero=settozero, package=package, cr=cr)$risk))
    print(i)
}

# Calculate ROC curves:

for (i in 1:7) {
    probs <- as.numeric(risks[,paste('model', i, sep='')])
    cutoffs <- sort(unique(probs), decreasing=TRUE)
    truepos <- rep(NA, length(cutoffs))
    falsepos <- rep(NA, length(cutoffs))
    auc <- rep(0.0, length(cutoffs))
    for (j in 1:length(cutoffs)) {
        ind <- as.numeric(probs > cutoffs[j])    
        truepos[j] <- sum(ind * agecir[,i])/sum(agecir[,i])
        falsepos[j] <- sum(ind * (1.0 - agecir[,i]))/sum(1.0 - agecir[,i])
        if (j > 1)
            auc[j] = (truepos[j] + truepos[j-1]) * (falsepos[j] - falsepos[j-1])
    }
    auc <- cumsum(auc) * 0.5
    roc <- cbind(cutoffs, truepos, falsepos, auc)
    save(roc, file=paste('ageroc', i, sep=''))
}

for (i in 1:7) {
    probs <- as.numeric(risks[,paste('model', i, sep='')])
    cutoffs <- sort(unique(probs), decreasing=TRUE)
    truepos <- rep(NA, length(cutoffs))
    falsepos <- rep(NA, length(cutoffs))
    auc <- rep(0.0, length(cutoffs))
    for (j in 1:length(cutoffs)) {
        ind <- as.numeric(probs > cutoffs[j])    
        truepos[j] <- sum(ind * cir[,i])/sum(cir[,i])
        falsepos[j] <- sum(ind * (1.0 - cir[,i]))/sum(1.0 - cir[,i])
        if (j > 1)
            auc[j] = (truepos[j] + truepos[j-1]) * (falsepos[j] - falsepos[j-1])
    }
    auc <- cumsum(auc) * 0.5
    roc <- cbind(cutoffs, truepos, falsepos, auc)
    save(roc, file=paste('roc', i, sep=''))
}

# Plot ROC curves:

# postscript(file.path(getwd(), 'rocs.eps'), paper='special', width=10, height=5, 
#            horizontal=FALSE)
op <- par(cex=1, mar=c(3.75,3.75,0.25,0.25), mfrow=c(1,2), mgp=c(2.5,1,0))

plot(1, xlim=c(0,1), ylim=c(0,1), type='n', xlab='False positive fraction', 
     ylab='True positive fraction')
abline(0, 1, lty='dashed')
cols=c('darkgray','red','blue','darkgreen','orange','purple','magenta')
aucs <- NULL
for (i in 1:7) {
    load(file=paste('roc', i, sep=''))
    aucs <- c(aucs, max(roc[,4]))
    lines(roc[,3], roc[,2], type='s', lwd=2, col=cols[i])
    for (j in c(0.05,0.1,0.15,0.2)) {
        tp <- approx(roc[,1], roc[,2], xout=j)$y
        fp <- approx(roc[,1], roc[,3], xout=j)$y
        idx <- nobs - findInterval(j,sort(roc[,1]))
            points(fp, tp, col=cols[i], pch=20)
        if (i == 1)
            text(fp, tp-0.015, labels=j, pos=4, offset=0.25, col=cols[i], 
                 cex=0.9)
    }
}
legend('bottomright', legend=paste('Model ', 1:7, '; AUC=',     
       format(round(aucs, 3), nsmall=3, scientific=FALSE), sep=''), 
        col=cols, lty=rep('solid',7), lwd=rep(2,7))

plot(1, xlim=c(0,1), ylim=c(0,1), type='n', xlab='False positive fraction', 
     ylab='True positive fraction')
abline(0, 1, lty='dashed')
cols=c('darkgray','red','blue','darkgreen','orange','purple','magenta')
aucs <- NULL
for (i in 1:7) {
    load(file=paste('ageroc', i, sep=''))
    aucs <- c(aucs, max(roc[,4]))
    lines(roc[,3], roc[,2], type='s', lwd=2, col=cols[i])
    for (j in c(0.05,0.1,0.15,0.2)) {
        tp <- approx(roc[,1], roc[,2], xout=j)$y
        fp <- approx(roc[,1], roc[,3], xout=j)$y
        idx <- nobs - findInterval(j,sort(roc[,1]))
        points(fp, tp, col=cols[i], pch=20)
        if (i == 1)
            text(fp, tp-0.015, labels=j, pos=4, offset=0.25, col=cols[i], 
                 cex=0.9)
    }
}
legend('bottomright', legend=paste('Model ', 1:7, '; AUC=',     
       format(round(aucs, 3), nsmall=3, scientific=FALSE), sep=''), 
        col=cols, lty=rep('solid',7), lwd=rep(2,7))

par(op)
# dev.off()
}
}
\keyword{datasets}
