/* ************************************************************************** */
/* File:            sampler.c                                                 */
/* Author:          Olli Saarela (olli.saarela@utoronto.ca                    */
/* Description:     Gibbs sampler and related functions.                      */
/* ************************************************************************** */

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>

#include "propadjust.h"
#include "mylib.h"

#define MAXTOTAL 10000

/* ************************************************************************** */

/* Global variables: */

int NITER, BURNIN, ADAPT, REFRESH, THIN, NOBS, NCOV, NCOVC, NAXS, BIRTHDEATH, NPPS, NPKG, TIMEVAR, SEED, COMP, NOTIME, LOGIT, NODATA;
double RHOA, RHOB, YEARS, DELTAI, DRANGE, sigma;

int steptotal;
int *steps, *maxstep, *pp, *predict, *include, *d, *pps, *pkg, *cr;
int **settozero;
long int **betacnt, **betaccnt, **sigmacnt;

double *sprob, *offset, *b, *smin, *smax, *vol, *y, *delta, *deltaold, *delta0, *delta0old, *deltamin, *deltamax,
*deltaminold, *deltamaxold, *rho, *beta, *betac;
double **z, **zc, **x, **bz, **bzold, **spos, **pred, **betasd, **betacsd, **sigmasd;
double ***lambda, ***lambdaold;

size_t ranks[MAXTOTAL];

gsl_rng *rgen;
gsl_ran_discrete_t *discd;

/* ************************************************************************** */

double loglik() {
    int i, j;
    double l=0.0;
    double mu[2];
    if (NODATA)
        return l;    
    if (LOGIT) {
        for (i = 0; i < NOBS; i++) {
            if (include[i]) {
                mu[0] = bz[0][i] + offset[i];
                mu[1] = bz[1][i];
                for (j = 0; j < NPKG; j++)
                    mu[cr[j]] += (*lambda[j][i]);
                l += (d[i] == 1) * mu[0] + (d[i] == 2) * mu[1] + (d[i] == 0) * log(sprob[i])
                     - log(sprob[i] + exp(mu[0]) + exp(mu[1]));
            }
        }
    }
    else {
        for (i = 0; i < NOBS; i++) {
            if (include[i]) {
                mu[0] = bz[0][i] + offset[i];
                mu[1] = bz[1][i];
                for (j = 0; j < NPKG; j++)
                    mu[cr[j]] += (*lambda[j][i]);
                l -= (y[i] - mu[0])*(y[i] - mu[0]);
            }            
        }
        l /= (2.0*sigma*sigma);
        l -= NOBS * log(sigma);
    }
    return l;
}

/* ************************************************************************** */

void update_bounds(int package) {
    int i;
    double u;
    double *dmin, *dmax;
    dmin = &(delta0[package]);
    dmax = &(delta0[package]);
    for (i = 0; i < steptotal; i++) {
        if (pkg[pp[i]] == package) {
            if (delta[i] > *dmax)
                dmax = &(delta[i]);
        }
    }
    u = gsl_ran_flat(rgen, *dmax - deltamax[package], *dmin - deltamin[package]);
    deltamin[package] += u;
    deltamax[package] += u;
    return;
}

/* ************************************************************************** */

void rescale(int package) {
    int i, counter=0;
    double mean=0.0;

    for (i = 0; i < NOBS; i++) {
        if (include[i]) {
            mean += (*lambda[package][i]);
            counter++;
        }
    }
    mean /= counter;
    delta0[package] -= mean;
    deltamin[package] -= mean;
    deltamax[package] -= mean;
    for (i = 0; i < steptotal; i++) {
        if (pkg[pp[i]] == package) {
            delta[i] -= mean;
        }
    }
    return;
}

/* ************************************************************************** */

double sumtozero(int package) {
    int i, counter=0;
    double mean=0.0;

    for (i = 0; i < NOBS; i++) {
        if (include[i]) {
            mean += (*lambda[package][i]);
            counter++;
        }
    }
    mean /= counter;
    return mean;
}

/* ************************************************************************** */

void savestate() {
    memcpy(deltaminold, deltamin, NPKG * sizeof(double));    
    memcpy(deltamaxold, deltamax, NPKG * sizeof(double));
    memcpy(deltaold, delta, steptotal * sizeof(double));
    memcpy(delta0old, delta0, NPKG * sizeof(double));
}

/* ************************************************************************** */

void restorestate() {
    memcpy(deltamin, deltaminold, NPKG * sizeof(double));    
    memcpy(deltamax, deltamaxold, NPKG * sizeof(double));
    memcpy(delta, deltaold, steptotal * sizeof(double));
    memcpy(delta0, delta0old, NPKG * sizeof(double));
}

/* ************************************************************************** */

int lowercorner(int id, int point) {
    int j;
    for (j = 0; j < NAXS; j++)
        if (spos[j][point] > x[j][id])
            break;
    return (j == NAXS);
}

/* ************************************************************************** */

void update_rho(int pid) {
    int i, axes=0;
    for (i = 0; i < NAXS; i++)
        axes += !(settozero[pid][i]);            
    rho[pid] = gsl_ran_gamma(rgen, RHOA + steps[pid], 1.0/(RHOB + vol[pid]));
    return;
}

/* ************************************************************************** */

int update_sigma(double sd) {

    /* Update sigma (random walk Metropolis): */

    double store, lold=loglik();

    store = sigma;
    sigma += gsl_ran_gaussian(rgen, sd);

    if (sigma < 0.0) {
        sigma = store;
        return 0;
    }
    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, exp(loglik() - lold))) {
        return 1;
    }
    else {
        sigma = store;
        return 0;
    }
}

/* ************************************************************************** */

int update_delta(int j) {

    /* Update delta (random walk/prior proposal Metropolis): */

    int i, k, package=pkg[pp[j]];
    double lold=loglik(), mhratio, store;
    double *dmin, *dmax;

    dmin = &(delta0[pkg[pp[j]]]);
    dmax = &(deltamax[pkg[pp[j]]]);
    for (i = 0; i < steptotal; i++) {
        if (pkg[pp[i]] == package) {
            if (i != j) {
                for (k = 0; k < NAXS; k++)
                    if (spos[k][i] < spos[k][j])
                        break;
                if (k == NAXS)
                    if (delta[i] < *dmax)
                        dmax = &(delta[i]);
                for (k = 0; k < NAXS; k++)
                    if (spos[k][i] > spos[k][j])
                        break;
                if (k == NAXS)
                    if (delta[i] > *dmin)
                        dmin = &(delta[i]);
            }
        }
    }
    savestate();
    store = delta[j];    
    memcpy(lambdaold[package], lambda[package], NOBS * sizeof(double *));
    delta[j] = gsl_ran_flat(rgen, fmax(delta[j] - DELTAI, *dmin), fmin(delta[j] + DELTAI, *dmax));

    if (delta[j] < store) {
        for (i = 0; i < NOBS; i++) {
            if (include[i]) {
                if (lambda[package][i] == &(delta[j])) {
                    lambda[package][i] = &(delta0[package]);
                    for (k = 0; k < steptotal; k++) {
                        if (pkg[pp[k]] == package)
                            if (delta[k] > *lambda[package][i])
                                if (lowercorner(i, k))
                                    lambda[package][i] = &(delta[k]);
                    }
                }
            }
        }
    }
    else if (delta[j] > store) {
        for (i = 0; i < NOBS; i++) {
            if (include[i]) {
                if (lambda[package][i] != &(delta[j])) {
                    if (delta[j] > *lambda[package][i])
                        if (lowercorner(i, j))
                            lambda[package][i] = &(delta[j]);
                }
            }
        }
    }
    rescale(package);
    mhratio = exp(loglik() - lold) *
              (gsl_ran_flat_pdf(deltaold[j], fmax(delta[j] - DELTAI, *dmin), fmin(delta[j] + DELTAI, *dmax)) /
               gsl_ran_flat_pdf(delta[j], fmax(deltaold[j] - DELTAI, *dmin), fmin(deltaold[j] + DELTAI, *dmax)));
    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, mhratio)) {
        return 1;
    }
    else {
        memcpy(lambda[package], lambdaold[package], NOBS * sizeof(double *));
        delta[j] = store;
        restorestate();
        return 0;
    }

}

/* ************************************************************************** */

int update_delta0(int package) {

    /* Update delta0 (random walk/prior proposal Metropolis): */

    int i;
    double lold=loglik(), mhratio, store;
    double *dmin, *dmax;

    dmin = &(deltamin[package]);
    dmax = &(deltamax[package]);
    for (i = 0; i < steptotal; i++) {
        if (pkg[pp[i]] == package)
            if (delta[i] < *dmax)
                dmax = &(delta[i]);
    }
    savestate();

    store = delta0[package];
    delta0[package] = gsl_ran_flat(rgen, fmax(delta0[package] - DELTAI, *dmin), fmin(delta0[package] + DELTAI, *dmax));
    rescale(package);

    mhratio = exp(loglik() - lold) *
              (gsl_ran_flat_pdf(store, fmax(delta0[package] - DELTAI, *dmin), fmin(delta0[package] + DELTAI, *dmax)) /
               gsl_ran_flat_pdf(delta0[package], fmax(store - DELTAI, *dmin), fmin(store + DELTAI, *dmax)));

    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, mhratio)) {
        return 1;
    }
    else {
        delta0[package] = store;
        restorestate();
        return 0;
    }
}

/* ************************************************************************** */

int update_beta(int j, double sd) {

    /* Update beta (random walk Metropolis): */

    int i;
    double store, lold=loglik(), mhratio;

    store = beta[j];
    beta[j] += gsl_ran_gaussian(rgen, sd);

    memcpy(bzold[0], bz[0], NOBS * sizeof(double));
    for (i = 0; i < NOBS; i++)
        bz[0][i] += (beta[j] - store) * z[j][i];

    mhratio = exp(loglik() - lold);
    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, mhratio)) {
        return 1;
    }
    else {
        beta[j] = store;
        memcpy(bz[0], bzold[0], NOBS * sizeof(double));
        return 0;
    }
}

/* ************************************************************************** */

int update_betac(int j, double sd) {

    /* Update betac (random walk Metropolis): */

    int i;
    double store, lold=loglik(), mhratio;

    store = betac[j];
    betac[j] += gsl_ran_gaussian(rgen, sd);

    memcpy(bzold[1], bz[1], NOBS * sizeof(double));
    for (i = 0; i < NOBS; i++)
        bz[1][i] += (betac[j] - store) * zc[j][i];

    mhratio = exp(loglik() - lold);
    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, mhratio)) {
        return 1;
    }
    else {
        betac[j] = store;
        memcpy(bz[1], bzold[1], NOBS * sizeof(double));

        return 0;
    }
}

/* ************************************************************************** */

int birth(int pid) {
    /* Propose new changepoint: */
    int i, j, package=pkg[pid];
    double lold=loglik(), mhratio;
    double *dmin, *dmax;

    pp[steptotal] = pid;

    for (i = 0; i < NAXS; i++) {
        if (settozero[pid][i])
            spos[i][steptotal] = 0.0;
        else {
            spos[i][steptotal] = gsl_ran_flat(rgen, smin[i], smax[i]);
        }
    }
    dmin = &(delta0[package]);
    dmax = &(deltamax[package]);
    for (i = 0; i < steptotal; i++) {
        if (pkg[pp[i]] == package) {
            for (j = 0; j < NAXS; j++)
                if (spos[j][i] < spos[j][steptotal])
                    break;            
            if (j == NAXS) 
                if (delta[i] < *dmax)
                    dmax = &(delta[i]);
            for (j = 0; j < NAXS; j++)
                if (spos[j][i] > spos[j][steptotal])
                    break;            
            if (j == NAXS)
                if (delta[i] > *dmin)
                    dmin = &(delta[i]);
        }
    }
    savestate();
    memcpy(lambdaold[package], lambda[package], NOBS * sizeof(double *));
    delta[steptotal] = gsl_ran_flat(rgen, *dmin, *dmax);

    for (i = 0; i < NOBS; i++) {
        if (include[i])
            if (delta[steptotal] > *lambda[package][i])
                if (lowercorner(i, steptotal))
                    lambda[package][i] = &(delta[steptotal]);
    }
    steptotal++;
    steps[pid]++;
    rescale(package);

    /* Accept/reject: */
    mhratio = exp(loglik() - lold) * rho[pid] * vol[pid] / steps[pid];
    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, mhratio)) {
        return 1;
    }
    else {
        memcpy(lambda[package], lambdaold[package], NOBS * sizeof(double *));
        steptotal--;
        steps[pid]--;
        restorestate();
        return 0;
    }
}

/* ************************************************************************** */

int death(int pid) {
    /* Propose to remove an existing changepoint: */
    int i=0, j=gsl_rng_uniform_int(rgen, steps[pid]), k=-1, l, package=pkg[pid];
    double lold=loglik(), mhratio, store;
    double *sold;
    sold = dvector(NAXS);

    while (i <= j) {
        k++;
        if (pp[k] == pid)
            i++;
    }
    j = k;

    savestate();
    store = delta[j];
    for (i = 0; i < NAXS; i++) {
        sold[i] = spos[i][j];
    }
    for (i = 0; i < NPKG; i++) {
        memcpy(lambdaold[i], lambda[i], NOBS * sizeof(double *));
    }
    for (i = 0; i < NOBS; i++) {
        if (include[i]) {
            if (lambda[package][i] == &(delta[j])) {
                lambda[package][i] = &(delta0[package]);
                for (k = 0; k < steptotal; k++) {
                    if (k != j) {
                        if (pkg[pp[k]] == package)
                            if (delta[k] > *lambda[package][i])
                                if (lowercorner(i, k))
                                    lambda[package][i] = &(delta[k]);
                    }
                }
            }
            for (k = 0; k < NPKG; k++) {
                if (lambda[k][i] != &(delta0[k])) {
                    if (lambda[k][i] > &(delta[j]))
                        lambda[k][i]--;
                }
            }
        }
    }
    steptotal--;
    steps[pid]--;
    for (i = j; i < steptotal; i++) {
        for (l = 0; l < NAXS; l++)
            spos[l][i] = spos[l][i+1];
        delta[i] = delta[i+1];
        pp[i] = pp[i+1];
    }
    rescale(package);

    /* Accept/reject: */
    mhratio = exp(loglik() - lold) / rho[pid] * (steps[pid] + 1) / vol[pid];
    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, mhratio)) {
        scrapdvector(sold);
        return 1;
    }
    else {
        for (i = 0; i < NPKG; i++) {
            memcpy(lambda[i], lambdaold[i], NOBS * sizeof(double *));
        }
        for (i = steptotal; i > j; i--) {
            for (l = 0; l < NAXS; l++)
                spos[l][i] = spos[l][i-1];
            delta[i] = delta[i-1];
            pp[i] = pp[i-1];
        }
        for (l = 0; l < NAXS; l++)
            spos[l][j] = sold[l];
        delta[j] = store;
        pp[j] = pid;
        steptotal++;
        steps[pid]++;
        restorestate();
        scrapdvector(sold);
        return 0;
    }
}

/* ************************************************************************** */

int move(int pid) {
    /* Propose to move an existing changepoint: */
    int i=0, j = gsl_rng_uniform_int(rgen, steps[pid]), k=-1, l, package=pkg[pid];

    while (i <= j) {
        k++;
        if (pp[k] == pid)
            i++;
    }
    j = k;

    double lold=loglik(), mhratio;
    double *sold;
    double **lowerlim, **upperlim;
    sold = dvector(NAXS);
    lowerlim = pdvector(NAXS);
    upperlim = pdvector(NAXS);

    for (i = 0; i < NAXS; i++) {
        sold[i] = spos[i][j];
        lowerlim[i] = &(smin[i]);
        upperlim[i] = &(smax[i]);
    }
    for (i = 0; i < steptotal; i++) {
        if (pkg[pp[i]] == package) {
            if (i != j) {
                for (l = 0; l < NAXS; l++) {
                    if (spos[l][i] <= spos[l][j] && spos[l][i] > *lowerlim[l])
                        lowerlim[l] = &(spos[l][i]);
                    if (spos[l][i] >= spos[l][j] && spos[l][i] < *upperlim[l])
                        upperlim[l] = &(spos[l][i]);
                }
            }
        }
    }
    for (i = 0; i < NAXS; i++) {
        if (settozero[pid][i])
            spos[i][j] = 0.0;
        else
            spos[i][j] = gsl_ran_flat(rgen, *lowerlim[i], *upperlim[i]);
    }

    savestate();
    memcpy(lambdaold[package], lambda[package], NOBS * sizeof(double *));
    for (i = 0; i < NOBS; i++) {
        if (include[i]) {
            if (lambda[package][i] == &(delta[j])) {
                lambda[package][i] = &(delta0[package]);
                for (k = 0; k < steptotal; k++) {
                    if (pkg[pp[k]] == package)
                        if (delta[k] > *lambda[package][i])
                            if (lowercorner(i, k))
                                lambda[package][i] = &(delta[k]);
                }
            }
            else if (lambda[package][i] != &(delta[j])) {
                if (delta[j] > *lambda[package][i])
                    if (lowercorner(i, j))
                        lambda[package][i] = &(delta[j]);
            }
        }
    }
    rescale(package);

    /* Accept/reject: */
    mhratio = exp(loglik() - lold);
    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, mhratio)) {
        scrapdvector(sold);
        scrappdvector(lowerlim);
        scrappdvector(upperlim);
        return 1;
    }
    else {
        memcpy(lambda[package], lambdaold[package], NOBS * sizeof(double *));
        for (i = 0; i < NAXS; i++)
            spos[i][j] = sold[i];
        restorestate();
        scrapdvector(sold);
        scrappdvector(lowerlim);
        scrappdvector(upperlim);
        return 0;
    }
}

/* ************************************************************************** */

int death_birth(int pid) {
    /* Combined death-birth move: */
    int i=0, j=gsl_rng_uniform_int(rgen, steps[pid]), k=-1, newpid, l, 
        package=pkg[pid], newpackage;
    double lold=loglik(), mhratio, store;
    double *min, *max;
    double *sold;
    sold = dvector(NAXS);

    while (i <= j) {
        k++;
        if (pp[k] == pid)
            i++;
    }
    j = k;
    newpid = pps[gsl_rng_uniform_int(rgen, NPPS)];
    /* newpid = gsl_ran_discrete(rgen, discd); */
    newpackage = pkg[newpid];

    /* Store old values: */

    savestate();
    store = delta[j];
    for (i = 0; i < NAXS; i++) {
        sold[i] = spos[i][j];
    }
    memcpy(lambdaold[package], lambda[package], NOBS * sizeof(double *));
    if (newpackage != package)
        memcpy(lambdaold[newpackage], lambda[newpackage], NOBS * sizeof(double *));

    /* Remove pointers to the removed point: */

    for (i = 0; i < NOBS; i++) {
        if (include[i]) {
            if (lambda[package][i] == &(delta[j])) {
                lambda[package][i] = &(delta0[package]);
                for (k = 0; k < steptotal; k++) {
                    if (k != j) {
                        if (pkg[pp[k]] == package)
                            if (delta[k] > *lambda[package][i])
                                if (lowercorner(i, k))
                                    lambda[package][i] = &(delta[k]);
                    }
                }
            }
        }
    }
    steps[pid]--;

    /* New point: */

    pp[j] = newpid;
    for (i = 0; i < NAXS; i++) {
        if (settozero[newpid][i])
            spos[i][j] = 0.0;
        else {
            spos[i][j] = gsl_ran_flat(rgen, smin[i], smax[i]);
        }
    }
    min = &(delta0[newpackage]);
    max = &(deltamax[newpackage]);
    for (i = 0; i < steptotal; i++) {
        if (pkg[pp[i]] == newpackage) {
            if (i != j) {
                for (l = 0; l < NAXS; l++)
                    if (spos[l][i] < spos[l][j])
                        break;
                if (l == NAXS)
                    if (delta[i] < *max)
                        max = &(delta[i]);
                for (l = 0; l < NAXS; l++)
                    if (spos[l][i] > spos[l][j])
                        break;
                if (l == NAXS)
                    if (delta[i] > *min)
                        min = &(delta[i]);
            }
        }
    }
    delta[j] = gsl_ran_flat(rgen, *min, *max);
    for (i = 0; i < NOBS; i++) {
        if (include[i]) {
            if (delta[j] > *lambda[newpackage][i])
                if (lowercorner(i, j))
                    lambda[newpackage][i] = &(delta[j]);
        }
    }
    steps[newpid]++;
    rescale(package);
    if (newpackage != package)
        rescale(newpackage);

    /* Accept/reject: */
    mhratio = exp(loglik() - lold) * (rho[newpid] / rho[pid]) *
              ((double)(steps[pid] + 1) / steps[newpid]) * (vol[newpid] / vol[pid]);
    if (gsl_rng_uniform_pos(rgen) < fmin(1.0, mhratio)) {
        scrapdvector(sold);
        return 1;
    }
    else {
        memcpy(lambda[package], lambdaold[package], NOBS * sizeof(double *));
        if (newpackage != package)
            memcpy(lambda[newpackage], lambdaold[newpackage], NOBS * sizeof(double *));
        steps[newpid]--;
        steps[pid]++;
        pp[j] = pid;
        for (i = 0; i < NAXS; i++)
            spos[i][j] = sold[i];
        delta[j] = store;
        restorestate();
        scrapdvector(sold);
        return 0;
    }
}

/* ************************************************************************** */

void getpred() {
    int i, j, k, l;
    double tlower, tupper, surv, bsurv, predsum=0.0;
    double hr[2], mu[2], ch[2];
    double **level;
    level = pdvector(NPKG);

    if (LOGIT) {
        if (NOTIME) {
            for (i = 0; i < NOBS; i++) {
                if (predict[i]) {
                    mu[0] = bz[0][i] + offset[i];
                    mu[1] = bz[1][i];
                    for (j = 0; j < NPKG; j++)
                        mu[cr[j]] += (*lambda[j][i]);
                    pred[i][0] = exp(mu[0])/(1.0 + exp(mu[0]) + exp(mu[1]));
                    pred[i][1] = exp(mu[1])/(1.0 + exp(mu[0]) + exp(mu[1]));
                }
            }
        }
        else {
            gsl_sort_index(ranks, spos[TIMEVAR], 1, steptotal);
            for (i = 0; i < NOBS; i++) {
                if (predict[i]) {
                    ch[0] = 0.0;
                    ch[1] = 0.0;
                    tlower = 0.0;
                    tupper = b[i];
                    pred[i][0] = 0.0;
                    pred[i][1] = 0.0;
                    surv = 1.0;
                    for (j = 0; j < NPKG; j++)
                        level[j] = &(delta0[j]);
                    for (j = 0; j < steptotal; j++) {
                        if (lowercorner(i, ranks[j])) {
                            if (spos[TIMEVAR][ranks[j]] < tupper) {
                                if (delta[ranks[j]] > *level[pkg[pp[ranks[j]]]]) {
                                    mu[0] = 0.0;
                                    mu[1] = 0.0;
                                    for (k = 0; k < NPKG; k++)
                                        mu[cr[k]] += (*level[k]);
                                    hr[0] = (exp(mu[0] + bz[0][i] + offset[i]))/(exp(mu[0] + bz[0][i] + offset[i]) + exp(mu[1] + bz[1][i]));
                                    hr[1] = (exp(mu[1] + bz[1][i]))/(exp(mu[0] + bz[0][i] + offset[i]) + exp(mu[1] + bz[1][i]));
                                    ch[0] += (spos[TIMEVAR][ranks[j]] - tlower) * exp(mu[0]);
                                    ch[1] += (spos[TIMEVAR][ranks[j]] - tlower) * exp(mu[1]);
                                    surv = exp(-(ch[0] * exp(bz[0][i] + offset[i]) + ch[1] * exp(bz[1][i])));
                                    tlower = spos[TIMEVAR][ranks[j]];
                                    level[pkg[pp[ranks[j]]]] = &(delta[ranks[j]]);
                                }
                            }
                            else
                                goto baseline;
                        }
                    }
                    baseline:
                    l = j;
                    mu[0] = 0.0;
                    mu[1] = 0.0;
                    for (k = 0; k < NPKG; k++)
                        mu[cr[k]] += (*level[k]);
                    hr[0] = (exp(mu[0] + bz[0][i] + offset[i]))/(exp(mu[0] + bz[0][i] + offset[i]) + exp(mu[1] + bz[1][i]));
                    hr[1] = (exp(mu[1] + bz[1][i]))/(exp(mu[0] + bz[0][i] + offset[i]) + exp(mu[1] + bz[1][i]));
                    ch[0] += (tupper - tlower) * exp(mu[0]);
                    ch[1] += (tupper - tlower) * exp(mu[1]);
                    surv = exp(-(ch[0] * exp(bz[0][i] + offset[i]) + ch[1] * exp(bz[1][i])));
                    bsurv = surv;
                    tlower = tupper;
                    tupper = b[i] + YEARS;
                    for (j = l; j < steptotal; j++) {
                        if (lowercorner(i, ranks[j])) {
                            if (spos[TIMEVAR][ranks[j]] < tupper) {
                                if (delta[ranks[j]] > *level[pkg[pp[ranks[j]]]]) {
                                    mu[0] = 0.0;
                                    mu[1] = 0.0;
                                    for (k = 0; k < NPKG; k++)
                                        mu[cr[k]] += (*level[k]);
                                    hr[0] = (exp(mu[0] + bz[0][i] + offset[i]))/(exp(mu[0] + bz[0][i] + offset[i]) + exp(mu[1] + bz[1][i]));
                                    hr[1] = (exp(mu[1] + bz[1][i]))/(exp(mu[0] + bz[0][i] + offset[i]) + exp(mu[1] + bz[1][i]));
                                    pred[i][0] += hr[0] * surv;

                                    pred[i][1] += hr[1] * surv;
                                    ch[0] += (spos[TIMEVAR][ranks[j]] - tlower) * exp(mu[0]);
                                    ch[1] += (spos[TIMEVAR][ranks[j]] - tlower) * exp(mu[1]);
                                    surv = exp(-(ch[0] * exp(bz[0][i] + offset[i]) + ch[1] * exp(bz[1][i])));
                                    pred[i][0] -= hr[0] * surv;
                                    pred[i][1] -= hr[1] * surv;
                                    tlower = spos[TIMEVAR][ranks[j]];
                                    level[pkg[pp[ranks[j]]]] = &(delta[ranks[j]]);
                                }
                            }
                            else
                                goto end;
                        }
                    }
                    end:
                    mu[0] = 0.0;
                    mu[1] = 0.0;
                    for (k = 0; k < NPKG; k++)
                        mu[cr[k]] += (*level[k]);
                    hr[0] = (exp(mu[0] + bz[0][i] + offset[i]))/(exp(mu[0] + bz[0][i] + offset[i]) + exp(mu[1] + bz[1][i]));
                    hr[1] = (exp(mu[1] + bz[1][i]))/(exp(mu[0] + bz[0][i] + offset[i]) + exp(mu[1] + bz[1][i]));
                    pred[i][0] += hr[0] * surv;
                    pred[i][1] += hr[1] * surv;
                    ch[0] += (tupper - tlower) * exp(mu[0]);
                    ch[1] += (tupper - tlower) * exp(mu[1]);
                    surv = exp(-(ch[0] * exp(bz[0][i] + offset[i]) + ch[1] * exp(bz[1][i])));
                    pred[i][0] -= hr[0] * surv;
                    pred[i][1] -= hr[1] * surv;
                    predsum += pred[i][0]/bsurv;
                    pred[i][0] /= bsurv;
                    pred[i][1] /= bsurv;
                }
            }
            Rprintf("Risksum=%f\n", predsum);
        }
    }
    else {
        for (i = 0; i < NOBS; i++) {
            if (predict[i]) {
                mu[0] = bz[0][i] + offset[i];
                mu[1] = bz[1][i];
                for (j = 0; j < NPKG; j++)
                    mu[cr[j]] += (*lambda[j][i]);
                pred[i][0] = mu[0];
                pred[i][1] = mu[1];
            }
        }
    }
    scrappdvector(level);
    return;
}

/* ************************************************************************** */

void sampler(int *iargs, double *dargs, int *idata, double *ddata,
             int *isettozero, int *ipackage, int *icr, int *steptotalf, int *stepsf,
             double *rhof, double *likf, double *betaf, double *betacf, double *lambdaf, double *predf, double *cpredf, double *sigmaf) {

    /* Copy constant arguments: */

    NITER = iargs[0];
    BURNIN = iargs[1];
    ADAPT = iargs[2];
    REFRESH = iargs[3];
    THIN = iargs[4];
    NOBS = iargs[5];
    NCOV = iargs[6];
    NCOVC = iargs[7];
    NAXS = iargs[8];
    BIRTHDEATH = iargs[9];
    NPPS = iargs[10];
    NPKG = iargs[11];
    TIMEVAR = iargs[12];
    SEED = iargs[13];
    COMP = iargs[14];
    NOTIME = iargs[15];
    LOGIT = iargs[16];
    NODATA = iargs[17];

    RHOA = dargs[0];
    RHOB = dargs[1];
    YEARS = dargs[2];
    DELTAI = dargs[3];
    DRANGE = dargs[4];

    /* Local variables: */

    int i, j, k, l, nupdates, counter, scounter=0, nsave=(NITER - BURNIN)/THIN;
    long int *deltaprop, *deltaacc, *delta0prop, *delta0acc;
    long int **prop, **acc;
    double *size;
    double sum;

    deltaprop = livector(NPKG);
    deltaacc = livector(NPKG);
    delta0prop = livector(NPKG);
    delta0acc = livector(NPKG);
    prop = limatrix(NPPS, 4);
    acc = limatrix(NPPS, 4);
    size = dvector(NPPS);

    gsl_set_error_handler_off();

    /* Counters to zero: */

    for (k = 0; k < NPKG; k++) {
        deltaprop[k] = 0;
        deltaacc[k] = 0;
        delta0prop[k] = 0;
        delta0acc[k] = 0;
    }
    for (k = 0; k < NPPS; k++) {
        for (l = 0; l < 4; l++) {
            prop[k][l] = 0;
            acc[k][l] = 0;            
        }
    }

    /* Initialize global variables: */

    steptotal = 0;

    steps = ivector(NPPS);
    maxstep = ivector(NPPS);
    pp = ivector(MAXTOTAL);
    predict = ivector(NOBS);
    include = ivector(NOBS);
    d = ivector(NOBS);
    pps = ivector(NPPS);
    pkg = ivector(NPPS);
    cr = ivector(NPKG); 

    settozero = imatrix(NPPS, NAXS);

    betacnt = limatrix(NCOV, 2);
    betaccnt = limatrix(NCOVC, 2);
    sigmacnt = limatrix(1, 2);

    sprob = dvector(NOBS);
    offset = dvector(NOBS);
    b = dvector(NOBS);
    y = dvector(NOBS);
    smin = dvector(NAXS);
    smax = dvector(NAXS);
    vol = dvector(NPPS);
    delta = dvector(MAXTOTAL);
    deltaold = dvector(MAXTOTAL);
    delta0 = dvector(NPKG);
    delta0old = dvector(NPKG);
    deltamin = dvector(NPKG);
    deltamax = dvector(NPKG);
    deltaminold = dvector(NPKG);
    deltamaxold = dvector(NPKG);
    rho = dvector(NPPS);
    beta = dvector(NCOV);
    betac = dvector(NCOVC);

    z = dmatrix(NCOV, NOBS);
    zc = dmatrix(NCOVC, NOBS);
    x = dmatrix(NAXS, NOBS);
    bz = dmatrix(2, NOBS);
    bzold = dmatrix(2, NOBS);
    spos = dmatrix(NAXS, MAXTOTAL);
    pred = dmatrix(NOBS, 2);
    betasd = dmatrix(NCOV, 3);
    betacsd = dmatrix(NCOVC, 3);
    sigmasd = dmatrix(1, 3);

    lambda = pdmatrix(NPKG, NOBS);
    lambdaold = pdmatrix(NPKG, NOBS);

    rgen = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rgen, SEED);

    /* Read data: */

    for (i = 0; i < NOBS; i++) {
        predict[i] = idata[vidx(i, 0, NOBS)];
        include[i] = idata[vidx(i, 1, NOBS)];
        d[i] = idata[vidx(i, 2, NOBS)];

        sprob[i] = ddata[vidx(i, 0, NOBS)];
        offset[i] = ddata[vidx(i, 1, NOBS)];
        b[i] = ddata[vidx(i, 2, NOBS)];
        y[i] = ddata[vidx(i, 3, NOBS)];

        counter = 4;
        for (j = 0; j < NAXS; j++) {
            x[j][i] = ddata[vidx(i, counter, NOBS)];
            counter++;
        }
        for (j = 0; j < NCOV; j++) {
            z[j][i] = ddata[vidx(i, counter, NOBS)];
            counter++;
        }
        for (j = 0; j < NCOVC; j++) {
            zc[j][i] = ddata[vidx(i, counter, NOBS)];
            counter++;
        }
    }

    Rprintf("Point process configuration:\n");
    for (i = 0; i < NPPS; i++) {
        pps[i] = i;
        pkg[i] = ipackage[i];
        Rprintf("%d ", pkg[i] + 1);
        for (j = 0; j < NAXS; j++) {
            settozero[i][j] = isettozero[vidx(i, j, NPPS)];
            Rprintf("%d ", settozero[i][j]);
        }
        Rprintf("\n");
    }
    Rprintf("  ");
    for (i = 0; i < NPPS; i++) {
        if (i == TIMEVAR)
            Rprintf("T ");
        else
            Rprintf(" ");
    }
    Rprintf("\n");
    Rprintf("Packages:\n");
    for (i = 0; i < NPKG; i++) {
        Rprintf("%d ", i + 1);
    }
    Rprintf("\n");
    for (i = 0; i < NPKG; i++) {
        cr[i] = icr[i];
        Rprintf("%d ", cr[i]);
    }
    Rprintf("\n");

    /* Set initial values: */

    for (i = 0; i < NAXS; i++) {
        smin[i] = 0.0;
        smax[i] = 1.0;
    }
    steptotal = 0;
    for (i = 0; i < NPPS; i++) {
        vol[i] = 1.0;
        steps[i] = 0;
        rho[i] = 1.0;
        maxstep[i] = NOBS;
    }
    for (i = 0; i < NPKG; i++) {
        delta0[i] = 0.0;
        deltamin[i] = -(double)DRANGE/2.0;
        deltamax[i] = (double)DRANGE/2.0;
    }
    for (i = 0; i < NCOV; i++) {
        beta[i] = 0.0;
        setupproposal(betacnt[i], betasd[i], 0.1);
    }
    beta[0] = betaf[0];
    for (i = 0; i < NCOVC; i++) {
        betac[i] = (COMP ? 0.0 : -1000.0);
        setupproposal(betaccnt[i], betacsd[i], 0.1);
    }
    betac[0] = (COMP ? betacf[0] : -1000.0);
    for (i = 0; i < NOBS; i++) {
        bz[0][i] = beta[0];
        bz[1][i] = betac[0];
    }
    for (i = 0; i < NOBS; i++) {
        for (j = 0; j < NPKG; j++) {
            lambda[j][i] = &(delta0[j]);
        }
    }
    sigma = sigmaf[0];
    setupproposal(sigmacnt[0], sigmasd[0], 0.1);

    /* Gibbs sampler loop: */

    time_t t0, t1, t2;
    t0 = time(NULL);
    t1 = time(NULL);
    double tdiff;

    Rprintf("Starting Gibbs sampler:\n");

    for (i = 1; i <= NITER; i++) {

        /* Adjust proposals: */

        if ((i < ADAPT && i % 100 == 0) || i == ADAPT) {
            for (k = 0; k < NCOV; k++)
                adjustproposal(betacnt[k], betasd[k], i, ADAPT, 0.1);
            if (COMP) {
                for (k = 0; k < NCOVC; k++)
                    adjustproposal(betaccnt[k], betacsd[k], i, ADAPT, 0.1);
            }
            if (!LOGIT)
                adjustproposal(sigmacnt[0], sigmasd[0], i, ADAPT, 0.1);
        }

        /* Birth/death moves: */

        for (k = 0; k < NPPS; k++)
            size[k] = rho[k] + 1.0;     
        discd = gsl_ran_discrete_preproc(NPPS, size);

        for (k = 0; k < BIRTHDEATH; k++) {
            l = pps[gsl_rng_uniform_int(rgen, NPPS)];
            /* l = gsl_ran_discrete(rgen, discd); */
            if (steps[l] == 0) {
                if (gsl_rng_uniform(rgen) < 0.5) {
                    prop[l][0]++;
                    acc[l][0] += birth(l);
                }
            }
            else if (steps[l] == maxstep[l]) {
                if (gsl_rng_uniform(rgen) < 0.5) {
                    prop[l][1]++;
                    acc[l][1] += death(l);
                }
            }
            else {
                if (gsl_rng_uniform(rgen) < 0.5) {
                    prop[l][0]++;
                    acc[l][0] += birth(l);
                }
                else {
                    prop[l][1]++;
                    acc[l][1] += death(l);
                }
            }
        }

        /* Combined death/birth moves: */

        for (k = 0; k < BIRTHDEATH; k++) {
            l = pps[gsl_rng_uniform_int(rgen, NPPS)];
            /* l = gsl_ran_discrete(rgen, discd); */
            if (steps[l] > 0) {
                prop[l][2]++;
                acc[l][2] += death_birth(l);
            }
        }

        /* Position change moves: */

        nupdates = (steptotal < 100) ? steptotal : 100;
        for (k = 0; k < nupdates; k++) {
            l = pp[gsl_rng_uniform_int(rgen, steptotal)];
            prop[l][3]++;
            acc[l][3] += move(l);
        }

        /* Other parameters: */

        nupdates = (steptotal+NPKG < 100) ? (steptotal+NPKG) : 100;
        for (k = 0; k < nupdates; k++) {
            l = gsl_rng_uniform_int(rgen, steptotal+NPKG);
            if (l >= steptotal) {
                delta0acc[l-steptotal] += update_delta0(l-steptotal);
                delta0prop[l-steptotal]++;
            }
            else {
                deltaacc[pkg[pp[l]]] += update_delta(l);
                deltaprop[pkg[pp[l]]]++;
            }
        }

        for (k = 0; k < NPPS; k++) {
            l = pps[gsl_rng_uniform_int(rgen, NPPS)];
            update_rho(l);
        }

        for (k = 0; k < NCOV; k++) {
            betacnt[k][0] += update_beta(k, betasd[k][0]);
            betacnt[k][1]++;
        }
        if (COMP) {
            for (k = 0; k < NCOVC; k++) {
                betaccnt[k][0] += update_betac(k, betacsd[k][0]);
                betaccnt[k][1]++;
            }
        }
        for (k = 0; k < NPKG; k++)
            update_bounds(k);

        if (!LOGIT) {
            sigmacnt[0][0] += update_sigma(sigmasd[0][0]);
            sigmacnt[0][1]++;
        }

        /* Write the output: */

        if (i > BURNIN && i % THIN == 0) {
            for (k = 0; k < NPKG; k++) {
                lambdaf[vidx(scounter * NPKG + k, 0, nsave * NPKG)] = k;
                counter = 0;
                for (j = 0; j < NOBS; j++) {
                    if (predict[j]) {
                        lambdaf[vidx(scounter * NPKG + k, counter + 1, nsave * NPKG)] = *lambda[k][j];
                        counter++;
                    }
                }
            }
            for (k = 0; k < NPPS; k++) {
                rhof[vidx(scounter, k, nsave)] = rho[k];
                stepsf[vidx(scounter, k, nsave)] = steps[k];
            }

            for (k = 0; k < NCOV; k++)
                betaf[vidx(scounter, k, nsave)] = beta[k];
            if (COMP) {
                for (k = 0; k < NCOVC; k++)
                    betacf[vidx(scounter, k, nsave)] = betac[k];
            }

            steptotalf[scounter] = steptotal;
            likf[scounter] = loglik();

            getpred();
            counter = 0;
            for (j = 0; j < NOBS; j++) {
                if (predict[j]) {
                    predf[vidx(scounter, counter, nsave)] = pred[j][0];
                    cpredf[vidx(scounter, counter, nsave)] = pred[j][1];
                    counter++;
                }
            }
            if (!LOGIT)            
                sigmaf[scounter] = sigma;
            scounter++;   
        }

        /* Time interval: */

        if (i % REFRESH == 0) {
            for (k = 0; k < NPKG; k++) {
                Rprintf("Acceptance ratio for delta0[%d] (=%f) updates: %f\n", k, delta0[k], (double)delta0acc[k]/delta0prop[k]);
            }
            for (l = 0; l < NPKG; l++) {
                sum = 0.0;
                counter = 0;
                for (k = 0; k < NOBS; k++) {
                    if (include[k]) {
                        sum += *lambda[l][k];
                        counter++;
                    }
                }
                Rprintf("mean=%f (%f, %f)\n", sum/(double)counter, deltamin[l], deltamax[l]);
            }
            t2 = time(NULL);
            tdiff = difftime(t2, t1);
            Rprintf("%d iterations (%d s)\n", i, (int)tdiff);
            Rprintf("%d points (", steptotal);
            for (j = 0; j < NPKG; j++) {
                l = 0;
                for (k = 0; k < steptotal; k++)
                    if (pkg[pp[k]] == j)
                        l++;
                Rprintf("%d", l);
                if (j < NPKG - 1)
                    Rprintf(", ");                
            }
            Rprintf("), l=%f\n", loglik());
            for (k = 0; k < NCOV; k++)
                Rprintf("%f ", beta[k]);
            Rprintf("\n");
            if (COMP) {
                for (k = 0; k < NCOVC; k++)
                    Rprintf("%f ", betac[k]);
            }
            Rprintf("\n");
            t1 = t2;
        }
    }
    t2 = time(NULL);
    tdiff = difftime(t2, t0);
    Rprintf("%d updates took %ld seconds.\n", (i-1), (long int)tdiff);

    for (l = 0; l < NPPS; l++) {
        Rprintf("Acceptance ratio for birth moves (%d): %f\n", l, (double)acc[l][0]/prop[l][0]);
        Rprintf("Acceptance ratio for death moves (%d): %f\n", l, (double)acc[l][1]/prop[l][1]);
        Rprintf("Acceptance ratio for combined death-birth moves (%d): %f\n", l, (double)acc[l][2]/prop[l][2]);
        Rprintf("Acceptance ratio for position change moves (%d): %f\n\n", l, (double)acc[l][3]/prop[l][3]);
    }
    for (k = 0; k < NPKG; k++) {
        Rprintf("Acceptance ratio for delta0[%d] updates: %f\n", k, (double)delta0acc[k]/delta0prop[k]);
        Rprintf("Acceptance ratio for delta[%d] updates: %f\n", k, (double)deltaacc[k]/deltaprop[k]);
    }
    for (k = 0; k < NCOV; k++)
        Rprintf("Acceptance ratio for beta[%d] updates (sd): %f (%f)\n", k, (double)betacnt[k][0]/betacnt[k][1], betasd[k][0]);
    if (COMP) {
        for (k = 0; k < NCOVC; k++)
            Rprintf("Acceptance ratio for betac[%d] updates (sd): %f (%f)\n", k, (double)betaccnt[k][0]/betaccnt[k][1], betacsd[k][0]);
    }
    if (!LOGIT)
        Rprintf("Acceptance ratio for sigma updates (sd): %f (%f)\n", (double)sigmacnt[0][0]/sigmacnt[0][1], sigmasd[0][0]);
    

    /* Free memory: */

    scraplivector(deltaprop); 
    scraplivector(deltaacc); 
    scraplivector(delta0prop);
    scraplivector(delta0acc);
    scraplimatrix(prop, NPPS);
    scraplimatrix(acc, NPPS);
    scrapdvector(size);

    scrapivector(steps);
    scrapivector(maxstep);
    scrapivector(pp);
    scrapivector(predict);
    scrapivector(include);
    scrapivector(d);
    scrapivector(pps);
    scrapivector(pkg);
    scrapivector(cr);    

    scrapimatrix(settozero, NPPS);

    scraplimatrix(betacnt, NCOV);
    scraplimatrix(betaccnt, NCOVC);
    scraplimatrix(sigmacnt, 1);

    scrapdvector(sprob);
    scrapdvector(offset);
    scrapdvector(b);
    scrapdvector(y);
    scrapdvector(smin);
    scrapdvector(smax);
    scrapdvector(vol);
    scrapdvector(delta);
    scrapdvector(deltaold);
    scrapdvector(delta0);
    scrapdvector(delta0old);
    scrapdvector(deltamin);
    scrapdvector(deltamax);
    scrapdvector(deltaminold);
    scrapdvector(deltamaxold);
    scrapdvector(rho);
    scrapdvector(beta);
    scrapdvector(betac);

    scrapdmatrix(z, NCOV);
    scrapdmatrix(zc, NCOVC);
    scrapdmatrix(x, NAXS);
    scrapdmatrix(bz, 2);
    scrapdmatrix(bzold, 2);
    scrapdmatrix(spos, NAXS);
    scrapdmatrix(pred, NOBS);
    scrapdmatrix(betasd, NCOV);
    scrapdmatrix(betacsd, NCOVC);
    scrapdmatrix(sigmasd, 1);

    scrappdmatrix(lambda, NPKG);
    scrappdmatrix(lambdaold, NPKG);

    gsl_rng_free(rgen);

    return;
}

/* ************************************************************************** */

static const R_CMethodDef CEntries[] = {
    {"sampler", (DL_FUNC) &sampler, 17},
    {NULL, NULL, 0}
};

void R_init_monoreg(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

/* ************************************************************************** */

