getcmat <- function(p) {
    cmat <- matrix(0, 2^p - 1, p)
    k <- 0
    for (i in 1:p) {
        comb <- t(combn(1:p, i))
        for (j in 1:nrow(comb)) {
            k <- k + 1
            cmat[k,] <- as.numeric(!(1:p %in% comb[j,]))
        }
    }
    return(cmat)
}

monoreg <- function(niter=15000, burnin=5000, adapt=5000, refresh=10, thin=5, birthdeath=10, seed=1,
                     rhoa=0.1, rhob=0.1, deltai=0.1, drange=2.0, predict, include, response, 
                     offset=NULL, axes, covariates, settozero, package) {
    
    nobs <- nrow(as.matrix(axes))
    naxs <- ncol(as.matrix(axes))
    ncov <- ncol(as.matrix(covariates))

    settozero <- as.matrix(settozero)
    npps <- nrow(settozero)
    npkg <- max(package)
    nsave <- (niter - burnin)/thin
    
    # These are redundant for linear regression:
    
    comp <- 0
    ccovariates <- rep(1.0, nobs)
    notime <- 1
    cr <- rep(0, length(unique(package)))
    timevar <- 0
    ncovc <- 1
    tstart <- rep(0.0, nobs)
    sprob <- rep(1.0, nobs)
    years <- 1.0
    
    # Check arguments:
    
    if (!(length(niter) == 1 & niter %% 1 == 0))
        stop('Argument niter must be integer-valued.')
    if (!(length(burnin) == 1 & burnin %% 1 == 0))
        stop('Argument burnin must be integer-valued.')
    if (!(length(adapt) == 1 & adapt %% 1 == 0))
        stop('Argument adapt must be integer-valued.')
    if (!(length(refresh) == 1 & refresh %% 1 == 0))
        stop('Argument refresh must be integer-valued.')
    if (!(length(thin) == 1 & thin %% 1 == 0))
        stop('Argument thin must be integer-valued.')
    if (!(length(birthdeath) == 1 & birthdeath %% 1 == 0))
        stop('Argument birthdeath must be integer-valued.')
    if (!(length(timevar) == 1 & timevar %% 1 == 0))
        stop('Argument timevar must be integer-valued.')
    if (!(length(seed) == 1 & seed %% 1 == 0))
        stop('Argument seed must be integer-valued.')
    if (!(length(rhoa) == 1 & rhoa > 0.0))
        stop('Argument rhoa must be a positive real value.')
    if (!(length(rhob) == 1 & rhob > 0.0))
        stop('Argument rhob must be a positive real value.')
    if (!(length(deltai) == 1 & deltai > 0.0))
        stop('Argument deltai must be a positive real value.')
    if (!(length(drange) == 1 & drange > 0.0))
        stop('Argument drange must be a positive real value.')
    if (length(package) != nrow(settozero))
        stop('Package number must be positive and cannot exceed the number of point processes.')
    if (min(package) < 1 | max(package) > npps)
        stop('Package number must be positive and cannot exceed the number of point processes.')
    if (!identical(as.integer(sort(unique(package))), 1:npkg))
        stop('The packages must be numbered from one to the total number of packages.')
    if (!identical(sort(union(as.numeric(settozero), c(0,1))), c(0,1)))
        stop('Argument settozero must involve a matrix of ones and zeros.')
    if (naxs != ncol(settozero))
        stop('Arguments settozero and axes do not match.')
    if (is.null(offset))
        offset <- rep(0.0, nobs)
    if (!(nrow(as.matrix(predict)) == nrow(as.matrix(include)) & 
          nrow(as.matrix(predict)) == nrow(as.matrix(response)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(sprob)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(axes)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(tstart)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(covariates)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(ccovariates)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(offset))))
        stop('Lengths of the data arguments do not match.')
    if (!(1 %in% colMeans(as.matrix(covariates))))
        stop('Argument covariates must involve a vector of ones to specify an intercept term.')
    if (!(1 %in% colMeans(as.matrix(ccovariates))))
        stop('Argument ccovariates must involve a vector of ones to specify an intercept term.')
    if (!(timevar %in% 0:naxs))
        stop('Argument timevar must be either zero or a column index of the argument axes.')
    if (min(axes) < 0 | max(axes) > 1)
        stop('The variables in argument axes must be scaled into zero-one interval.')
    if (nsave < 1 | nsave > niter)
        stop('Number of interations to be saved is too small or too large.')
    if (nsave %% 1 != 0)
        stop('Number of interations to be saved is not integer-valued.')
    if (niter <= adapt)
        stop('Argument adapt must be smaller than argument niter.')    
    if (niter <= burnin)
        stop('Argument burnin must be smaller than argument niter.')
    if (adapt < 1000)
        warning('Use at least 1000 iterations for adapting the Metropolis-Hastings proposals.')
    if (!identical(sort(union(include, c(0,1))), c(0,1)))
        stop('Argument include must be a vector of ones and zeros.')
    if (sum(include) == 0)
        stop('Argument include must specify observations to be included in the likelihood.')
    if (!identical(sort(union(predict, c(0,1))), c(0,1)))
        stop('Argument predict must be a vector of ones and zeros.')
    if (sum(predict) == 0)
        stop('Argument predict must specify observations for prediction.')
    
    results <- .C("sampler", 
                  as.integer(c(niter, burnin, adapt, refresh, thin, nobs, ncov, ncovc, naxs, birthdeath, npps, npkg, timevar-1, seed, comp, notime, 0, 0)),
                  as.double(c(rhoa, rhob, years, deltai, drange)),
                  as.integer(cbind(predict, include, response)),
                  as.double(cbind(sprob, offset, tstart, response, axes, covariates, ccovariates)),
                  as.integer(settozero),
                  as.integer(package - 1),
                  as.integer(cr),
                  as.integer(numeric(nsave)),
                  as.integer(numeric(nsave * npps)),
                  numeric(nsave * npps),
                  numeric(nsave),
                  as.numeric(rep(mean(response), nsave * ncov)),
                  numeric(nsave * ncovc),
                  numeric(nsave * npkg * (sum(predict) + 1)),
                  numeric(nsave * sum(predict)),
                  numeric(nsave * sum(predict)),
                  as.numeric(rep(sd(response), nsave)),
                  PACKAGE='monoreg')
    
    return(list('steptotal'=results[[8]], 'steps'=matrix(results[[9]], nsave, npps), 'rho'=matrix(results[[10]], nsave, npps),
                'loglik'=results[[11]], 'beta'=matrix(results[[12]], nsave, ncov),
                'phi'=matrix(results[[14]], nsave * npkg, sum(predict) + 1), 'pred'=matrix(results[[15]], nsave, sum(predict)),
                'sigma'=results[[17]]))
}

monosurv <- function(niter=15000, burnin=5000, adapt=5000, refresh=10, thin=5, birthdeath=10, timevar=0, seed=1,
               rhoa=0.1, rhob=0.1, years=NULL, deltai=0.1, drange=2.0, predict, include, casestatus, sprob=NULL, 
               offset=NULL, tstart=NULL, axes, covariates, ccovariates=NULL, settozero, package, cr=NULL) {

    nobs <- nrow(as.matrix(axes))
    naxs <- ncol(as.matrix(axes))
    ncov <- ncol(as.matrix(covariates))
    {
    if (!is.null(ccovariates))
        ncovc <- ncol(as.matrix(ccovariates))
    else
        ncovc <- 1
    }
    settozero <- as.matrix(settozero)
    npps <- nrow(settozero)
    npkg <- max(package)
    nsave <- (niter - burnin)/thin
    comp <- 1
    notime <- ifelse(timevar == 0, 1, 0)
    if (is.null(cr))
        cr <- rep(0, length(unique(package)))

    # Check arguments:

    if (!(length(niter) == 1 & niter %% 1 == 0))
        stop('Argument niter must be integer-valued.')
    if (!(length(burnin) == 1 & burnin %% 1 == 0))
        stop('Argument burnin must be integer-valued.')
    if (!(length(adapt) == 1 & adapt %% 1 == 0))
        stop('Argument adapt must be integer-valued.')
    if (!(length(refresh) == 1 & refresh %% 1 == 0))
        stop('Argument refresh must be integer-valued.')
    if (!(length(thin) == 1 & thin %% 1 == 0))
        stop('Argument thin must be integer-valued.')
    if (!(length(birthdeath) == 1 & birthdeath %% 1 == 0))
        stop('Argument birthdeath must be integer-valued.')
    if (!(length(timevar) == 1 & timevar %% 1 == 0))
        stop('Argument timevar must be integer-valued.')
    if (!(length(seed) == 1 & seed %% 1 == 0))
        stop('Argument seed must be integer-valued.')
    if (!(length(rhoa) == 1 & rhoa > 0.0))
        stop('Argument rhoa must be a positive real value.')
    if (!(length(rhob) == 1 & rhob > 0.0))
        stop('Argument rhob must be a positive real value.')
    if (!(length(deltai) == 1 & deltai > 0.0))
        stop('Argument deltai must be a positive real value.')
    if (!(length(drange) == 1 & drange > 0.0))
        stop('Argument drange must be a positive real value.')
    if (length(package) != nrow(settozero))
        stop('Package number must be positive and cannot exceed the number of point processes.')
    if (min(package) < 1 | max(package) > npps)
        stop('Package number must be positive and cannot exceed the number of point processes.')
    if (!identical(as.integer(sort(unique(package))), 1:npkg))
        stop('The packages must be numbered from one to the total number of packages.')
    if (!identical(sort(union(as.numeric(settozero), c(0,1))), c(0,1)))
        stop('Argument settozero must involve a matrix of ones and zeros.')
    if (naxs != ncol(settozero))
        stop('Arguments settozero and axes do not match.')
    if (!identical(sort(union(cr, c(0,1))), c(0,1)))
        stop('Argument cr must be a vector of ones and zeros.')
    if (!identical(sort(union(casestatus, c(0,1,2))), c(0,1,2)))
        stop('Argument casestatus must only involve numbers zero, one or two')
    if (!(2 %in% casestatus) & !is.null(ccovariates))
        stop('Argument casestatus must involve competing events.')
    if (!(2 %in% casestatus) & max(cr) == 1)
        stop('Argument casestatus must involve competing events.')
    if (2 %in% casestatus & !(1 %in% casestatus))
        stop('Argument casestatus must involve events of interest.')
    if (2 %in% casestatus & is.null(ccovariates))
        stop('Argument ccovariates must be specified (to include at least an intercept term).')
    if (is.null(sprob))
        sprob <- rep(1.0, nobs)
    if (is.null(tstart))
        tstart <- rep(0.0, nobs)
    if (min(tstart) < 0 | max(tstart) > 1)
        stop('The variable in argument tstart must be scaled into zero-one interval.')
    if (is.null(offset))
        offset <- rep(0.0, nobs)
    if (!(2 %in% casestatus) & is.null(ccovariates)) {
        comp <- 0
        ccovariates <- rep(1.0, nobs)
    }
    if (!(nrow(as.matrix(predict)) == nrow(as.matrix(include)) & 
          nrow(as.matrix(predict)) == nrow(as.matrix(casestatus)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(sprob)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(axes)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(tstart)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(covariates)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(ccovariates)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(offset))))
        stop('Lengths of the data arguments do not match.')
    if (!(1 %in% colMeans(as.matrix(covariates))))
        stop('Argument covariates must involve a vector of ones to specify an intercept term.')
    if (!(1 %in% colMeans(as.matrix(ccovariates))))
        stop('Argument ccovariates must involve a vector of ones to specify an intercept term.')
    if (!(timevar %in% 0:naxs))
        stop('Argument timevar must be either zero or a column index of the argument axes.')
    if (timevar %in% 1:naxs & is.null(years))
        stop('Argument years must be specified.')    
    if (timevar == 0 & !is.null(years))
        stop('Argument years is specified, so argument timevar must specify a time variable.')
    if (is.null(years))
        years <- 1.0
    if (!(length(years) == 1 & years > 0.0 & years <= 1.0))
        stop('Argument years must be a real value in the interval (0,1].')
    if (min(axes) < 0 | max(axes) > 1)
        stop('The variables in argument axes must be scaled into zero-one interval.')
    if (nsave < 1 | nsave > niter)
        stop('Number of interations to be saved is too small or too large.')
    if (nsave %% 1 != 0)
        stop('Number of interations to be saved is not integer-valued.')
    if (niter <= burnin)
        stop('Argument burnin must be smaller than argument niter.')
    if (niter <= adapt)
        stop('Argument adapt must be smaller than argument niter.')    
    if (adapt < 1000)
        warning('Use at least 1000 iterations for adapting the Metropolis-Hastings proposals.')
    if (!identical(sort(union(include, c(0,1))), c(0,1)))
        stop('Argument include must be a vector of ones and zeros.')
    if (sum(include) == 0)
        stop('Argument include must specify observations to be included in the likelihood.')
    if (!identical(sort(union(predict, c(0,1))), c(0,1)))
        stop('Argument predict must be a vector of ones and zeros.')
    if (sum(predict) == 0)
        stop('Argument predict must specify new observations for prediction.')

    p1 <- mean(casestatus == 1)
    p2 <- mean(casestatus == 2)
    
    results <- .C("sampler", 
    as.integer(c(niter, burnin, adapt, refresh, thin, nobs, ncov, ncovc, naxs, birthdeath, npps, npkg, timevar-1, seed, comp, notime, 1, 0)),
    as.double(c(rhoa, rhob, years, deltai, drange)),
    as.integer(cbind(predict, include, casestatus)),
    as.double(cbind(sprob, offset, tstart, casestatus, axes, covariates, ccovariates)),
    as.integer(settozero),
    as.integer(package - 1),
    as.integer(cr),
    as.integer(numeric(nsave)),
    as.integer(numeric(nsave * npps)),
    numeric(nsave * npps),
    numeric(nsave),
    as.numeric(rep(ifelse(p1 > 0, log(1+p2/(1.0-p2))-log((1.0-p1)/p1 - p2/(1.0-p2)), 0.0), nsave * ncov)),
    as.numeric(rep(ifelse(p2 > 0, log(1+p1/(1.0-p1))-log((1.0-p2)/p2 - p1/(1.0-p1)), 0.0), nsave * ncovc)),
    numeric(nsave * npkg * (sum(predict) + 1)),
    numeric(nsave * sum(predict)),
    numeric(nsave * sum(predict)),
    as.numeric(rep(sd(casestatus), nsave)),
    PACKAGE='monoreg')

    return(list('steptotal'=results[[8]], 'steps'=matrix(results[[9]], nsave, npps), 'rho'=matrix(results[[10]], nsave, npps),
                'loglik'=results[[11]], 'beta'=matrix(results[[12]], nsave, ncov), 'betac'=matrix(results[[13]], nsave, ncovc),
                'phi'=matrix(results[[14]], nsave * npkg, sum(predict) + 1), 'risk'=matrix(results[[15]], nsave, sum(predict)),
                'crisk'=matrix(results[[16]], nsave, sum(predict))))
}

ordmonoreg <- function(niter=15000, burnin=5000, adapt=5000, refresh=10, thin=20, 
                       birthdeath=1, logit=FALSE, gam=FALSE, seed=1, rhoa=0.1, rhob=0.1, 
                       deltai=0.5, dlower=0, dupper=1, invprob=1.0, dc=0.0,
                       predict, include, outcome, axes, covariates=NULL, 
                       cluster=NULL, ncluster=NULL, settozero) {
    nobs <- nrow(as.matrix(axes))
    naxs <- ncol(as.matrix(axes))
    {
        if (!is.null(covariates)) {
            if (!is.null(cluster))
                covariates <- cbind(cluster, covariates)
            else
                covariates <- cbind(1, covariates)
            ncov <- ncol(as.matrix(covariates))
        }
        else {
            if (!is.null(cluster))
                covariates <- cluster
            else {
                covariates <- rep(1, nobs)
            }
            ncov <- 1
        }
        if (is.null(cluster))
            ncluster <- 1;    
    }
    settozero <- as.matrix(settozero)
    npps <- nrow(settozero)
    nsave <- (niter - burnin)/thin
    ncat <- length(unique(outcome))
    
    {if (!is.logical(logit))
        stop('Argument logit must be logical.')
        else
            logit <- as.integer(logit)
    }
    {if (logit == 1) {
        if (!(length(dlower) == 1 & deltai < 0.0))
            stop('Argument dlower must be a negative real value.')    
        if (!(length(deltai) == 1 & deltai > 0.0))
            stop('Argument dupper must be a positive real value.')    
    }
        else {
            if (!(length(dlower) == 1 & dlower >= 0.0 & dlower <= 1.0))
                stop('Argument dlower must be between zero and one.')    
            if (!(length(dupper) == 1 & dupper >= 0.0 & dupper <= 1.0))
                stop('Argument dupper must be between zero and one.')    
        }}
    if (dlower >= dupper)
        stop('Argument dupper must be greater than dlower.')
    if (!(length(invprob) == 1 & invprob >= 0.0 & invprob <= 1.0))
        stop('Argument invprob must be between zero and one.')    
    {if (!is.logical(gam))
        stop('Argument gam must be logical.')
        else
            gam <- as.integer(gam)
        }
    nint <- 1
    if (gam==1 & ncov > 1) {
        for (i in 2:ncov)
            nint <- ifelse(length(unique(covariates[,i])) > nint, length(unique(covariates[,i])), nint)
    }
    if (!(length(niter) == 1 & niter %% 1 == 0))
        stop('Argument niter must be integer-valued.')
    if (!(length(burnin) == 1 & burnin %% 1 == 0))
        stop('Argument burnin must be integer-valued.')
    if (!(length(adapt) == 1 & adapt %% 1 == 0))
        stop('Argument adapt must be integer-valued.')
    if (!(length(refresh) == 1 & refresh %% 1 == 0))
        stop('Argument refresh must be integer-valued.')
    if (!(length(thin) == 1 & thin %% 1 == 0))
        stop('Argument thin must be integer-valued.')
    if (!(length(birthdeath) == 1 & birthdeath %% 1 == 0))
        stop('Argument birthdeath must be integer-valued.')
    if (!(length(seed) == 1 & seed %% 1 == 0))
        stop('Argument seed must be integer-valued.')
    if (!(length(rhoa) == 1 & rhoa > 0.0))
        stop('Argument rhoa must be a positive real value.')
    if (!(length(rhob) == 1 & rhob > 0.0))
        stop('Argument rhob must be a positive real value.')
    if (!(length(deltai) == 1 & deltai > 0.0))
        stop('Argument deltai must be a positive real value.')
    if (!(length(dc) == 1 & dc >= 0.0))
        stop('Argument dc must be a non-negative real value.')      
    if (!identical(sort(union(as.integer(settozero), c(0,1))), c(0,1)))
        stop('Argument settozero must involve a matrix of ones and zeros.')
    if (naxs != ncol(settozero))
        stop('Arguments settozero and axes do not match.')
    if (!identical(sort(union(outcome, seq(1,ncat,1))), seq(1,ncat,1)))
        stop('Argument outcome must only involve integer numbers')
    if (!(nrow(as.matrix(predict)) == nrow(as.matrix(include)) & 
          nrow(as.matrix(predict)) == nrow(as.matrix(outcome)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(axes)) &
          nrow(as.matrix(predict)) == nrow(as.matrix(covariates))))
        stop('Lengths of the data arguments do not match.')
    if (min(axes) < 0 | max(axes) > 1)
        stop('The variables in argument axes must be scaled into zero-one interval.')
    if (nsave < 1 | nsave > niter)
        stop('Number of interations to be saved is too small or too large.')
    if (nsave %% 1 != 0)
        stop('Number of interations to be saved is not integer-valued.')
    if (niter <= burnin)
        stop('Argument burnin must be smaller than argument niter.')
    if (!identical(sort(union(include, c(0,1))), c(0,1)))
        stop('Argument include must be a vector of ones and zeros.')
    if (sum(include) == 0)
        stop('Argument include must specify observations to be included in the likelihood.')
    if (!identical(sort(union(predict, c(0,1))), c(0,1)))
        stop('Argument predict must be a vector of ones and zeros.')
    if (sum(predict) == 0)
        stop('Argument predict must specify new observations for prediction.')
    if (!is.null(cluster)) {
        if(sum(sort(unique(cluster)) == seq(1:ncluster) - 1) != ncluster)
            stop('Clusters must be numbered from 0 to ncluster-1.')
    }
    
    results <- .C("ordsampler", 
                  as.integer(c(niter, burnin, adapt, refresh, thin, nobs, ncov, naxs, birthdeath, npps, ncat, ncluster, logit, seed, gam)),
                  as.double(c(rhoa, rhob, deltai, dlower, dupper, invprob, dc)),
                  as.integer(cbind(predict, include, outcome-1)),
                  as.double(cbind(axes, covariates)),
                  as.integer(settozero),
                  as.integer(numeric(nsave)),
                  as.integer(numeric(nsave * npps)),
                  numeric(nsave * npps),
                  numeric(nsave),
                  numeric(nsave),
                  numeric(nsave * ncluster),
                  numeric(nsave * ncov * nint),
                  numeric(nsave * ncat * (sum(predict) + 1)),
                  numeric(nsave * ncov),
                  PACKAGE='monoreg')
    
    return(list('steptotal'=results[[6]], 'steps'=matrix(results[[7]], nsave, npps), 'rho'=matrix(results[[8]], nsave, npps),
                'loglik'=results[[9]], 'tau'=results[[10]], 'alpha'=matrix(results[[11]], nsave, ncluster),
                'beta'=array(results[[12]], c(nsave, ncov, nint)),
                'lambda'=matrix(results[[13]], nsave * ncat, sum(predict) + 1),
                'sigmasq'=matrix(results[[14]], nsave, ncov)))
}
