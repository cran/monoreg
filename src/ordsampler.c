/* ************************************************************************** */
/* File:            ordsampler.c                                              */
/* Author:          Olli Saarela (olli.saarela@utoronto.ca                    */
/* Description:     Gibbs sampler for ordinal monotonic regression.           */
/* ************************************************************************** */

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>

#include "propadjust.h"
#include "mylib.h"

#define MAXTOTAL 10000

/* ************************************************************************** */

/* Global variables: */

extern int NITER, BURNIN, ADAPT, REFRESH, THIN, NOBS, NCOV, NAXS, BIRTHDEATH, NPPS, SEED, LOGIT, steptotal;
int NCAT, NCLUSTER, GAM, ninc;

extern double RHOA, RHOB, DELTAI;
double DLOWER, DUPPER, INVPROB, tau, DC;

extern int *steps, *maxstep, *pp, *predict, *include, *d, *pps;
int *dimtotals, *direction, *nint;

extern int **settozero;
extern long int **betacnt;
long int **alphacnt;

extern double *smin, *smax, *vol, *delta0, *delta0old, *deltamin, *deltamax,
*deltaminold, *deltamaxold, *rho;

double *sigmasq, *bz1, *bz1old, *alpha;

extern double **z, **x, **spos, **betasd;
double **alphasd, **delta1, **delta1old, **betao;

extern double ***lambda, ***lambdaold;

extern size_t ranks[MAXTOTAL];

extern gsl_rng *rgen;
extern gsl_ran_discrete_t *discd;

/* ************************************************************************** */

double ordloglik() {
    int i;
    double l=0.0;
	if (LOGIT) {
		for (i = 0; i < NOBS; i++) {
			if (include[i]) {
				if (d[i] == NCAT-1)
					l += log(1.0/(1.0 + exp(-(*lambda[d[i]][i] + bz1[i]))));
				else if (d[i] == 0)
					l += log(1.0 - 1.0/(1.0 + exp(-(*lambda[d[i]+1][i] + bz1[i]))));
				else
					l += log(1.0/(1.0 + exp(-(*lambda[d[i]][i] + bz1[i]))) - 1.0/(1.0 + exp(-(*lambda[d[i]+1][i] + bz1[i]))));
			}
		}
	}
	else {
		for (i = 0; i < NOBS; i++) {
			if (include[i]) {
				if (d[i] == NCAT-1)
					l += log(*lambda[d[i]][i]);
				else if (d[i] == 0)
					l += log(1.0 - *lambda[d[i]+1][i]);	
				else
					l += log(*lambda[d[i]][i] - *lambda[d[i]+1][i]);
			}
		}
	}
    return l;
}

/* ************************************************************************** */

double ordloglikcat(int j) {
    int i;
    double l=0.0;
	if (LOGIT) {
		for (i = 0; i < NOBS; i++) {		
			if (include[i]) {
				if (z[0][i] == j) {
					if (d[i] == NCAT-1)
						l += log(1.0/(1.0 + exp(-(*lambda[d[i]][i] + bz1[i]))));
					else if (d[i] == 0)
						l += log(1.0 - 1.0/(1.0 + exp(-(*lambda[d[i]+1][i] + bz1[i]))));
					else
						l += log(1.0/(1.0 + exp(-(*lambda[d[i]][i] + bz1[i]))) - 1.0/(1.0 + exp(-(*lambda[d[i]+1][i] + bz1[i]))));
				}
			}
		}
	}
    return l;
}

/* ************************************************************************** */

void ordsavestate() {
    int level;
    memcpy(deltaminold, deltamin, NCAT * sizeof(double));    
    memcpy(deltamaxold, deltamax, NCAT * sizeof(double));
    memcpy(delta0old, delta0, NCAT * sizeof(double));
    for (level = 0; level < NCAT; level++)
        memcpy(delta1old[level], delta1[level], steptotal * sizeof(double));
}

/* ************************************************************************** */

void ordrestorestate() {
    int level;
    memcpy(deltamin, deltaminold, NCAT * sizeof(double));    
    memcpy(deltamax, deltamaxold, NCAT * sizeof(double));
    memcpy(delta0, delta0old, NCAT * sizeof(double));
    for (level = 0; level < NCAT; level++)
        memcpy(delta1[level], delta1old[level], steptotal * sizeof(double));
}

/* ************************************************************************** */

int ordlowercorner(int id, int point) {
    int j;
    for (j = 0; j < NAXS; j++)
        if (spos[j][point] > x[j][id])
            break;
    return (j == NAXS);
}

/* ************************************************************************** */

void ordinvert(int k) {
	int i;
	for (i = 0; i < NOBS; i++) {
		x[k][i] = 1.0 - x[k][i];
	}
	return;
}

/* ************************************************************************** */

int ordfindpoint(int id, int level) {
	int k, point=-1;
	double d=delta0[level];
	for (k = 0; k < steptotal; k++) {
		if (delta1[level][k] > d) {
			if (ordlowercorner(id, k)) {
				d = delta1[level][k];
				point=k;
			}
		}
	}
	return point;
}

/* ************************************************************************** */

void ordupdate_rho(int pid) {
    int i, axes=0;
    for (i = 0; i < NAXS; i++)
        axes += !(settozero[pid][i]);            
    rho[pid] = gsl_ran_gamma(rgen, RHOA + steps[pid], 1.0/(RHOB + vol[pid]));
    return;
}

/* ************************************************************************** */

int ordupdate_delta(int j, int level) {

    /* Update delta (random walk/prior proposal Metropolis): */

    int i, k;
    double lold=ordloglik(), mhratio, store;
    double *dmin, *dmax;

    dmin = &(delta0[level]);
    dmax = &(deltamax[level]);
    for (i = 0; i < steptotal; i++) {
        if (i != j) {
            for (k = 0; k < NAXS; k++)
                if (spos[k][i] < spos[k][j])
                    break;
            if (k == NAXS)
                if (delta1[level][i] < *dmax)
                    dmax = &(delta1[level][i]);
            for (k = 0; k < NAXS; k++)
                if (spos[k][i] > spos[k][j])
                    break;
            if (k == NAXS)
                if (delta1[level][i] > *dmin)
                    dmin = &(delta1[level][i]);
        }
    }
    if (delta1[level-1][j] < *dmax)
        dmax = &(delta1[level-1][j]);
    if (level < NCAT-1) {
        if (delta1[level+1][j] > *dmin)
            dmin = &(delta1[level+1][j]);
    }
    ordsavestate();
    store = delta1[level][j];
    memcpy(lambdaold[level], lambda[level], NOBS * sizeof(double *));
    delta1[level][j] = gsl_ran_flat(rgen, fmax(delta1[level][j] - DELTAI, *dmin), fmin(delta1[level][j] + DELTAI, *dmax));

    if (delta1[level][j] < store) {
        for (i = 0; i < NOBS; i++) {
			if (lambda[level][i] == &(delta1[level][j])) {
				lambda[level][i] = &(delta0[level]);
				for (k = 0; k < steptotal; k++) {
					if (delta1[level][k] > *lambda[level][i])
						if (ordlowercorner(i, k))
							lambda[level][i] = &(delta1[level][k]);
				}
			}
        }
    }
    else if (delta1[level][j] > store) {
        for (i = 0; i < NOBS; i++) {
			if (lambda[level][i] != &(delta1[level][j])) {
				if (delta1[level][j] > *lambda[level][i])
					if (ordlowercorner(i, j))
						lambda[level][i] = &(delta1[level][j]);
			}
        }
    }

    mhratio = exp(ordloglik() - lold) *
              (gsl_ran_flat_pdf(delta1old[level][j], fmax(delta1[level][j] - DELTAI, *dmin), fmin(delta1[level][j] + DELTAI, *dmax)) /
               gsl_ran_flat_pdf(delta1[level][j], fmax(delta1old[level][j] - DELTAI, *dmin), fmin(delta1old[level][j] + DELTAI, *dmax)));
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        return 1;
    }
    else {
        memcpy(lambda[level], lambdaold[level], NOBS * sizeof(double *));
        delta1[level][j] = store;
        ordrestorestate();
        return 0;
    }

}

/* ************************************************************************** */

int ordupdate_delta_joint(int j) {

    /* Update delta jointly (random walk/prior proposal Metropolis): */

    int i, k, level, order;
    double lold=ordloglik(), mhratio;
	double *store;	
    double **dmin, **dmax;

    store = dvector(NCAT);
    dmin = pdvector(NCAT);
    dmax = pdvector(NCAT);

    for (level = 1; level < NCAT; level++) {
        dmin[level] = &(delta0[level]);
        dmax[level] = &(deltamax[level]);
        for (i = 0; i < steptotal; i++) {
			if (i != j) {
				for (k = 0; k < NAXS; k++)
					if (spos[k][i] < spos[k][j])
						break;            
				if (k == NAXS) 
					if (delta1[level][i] < *dmax[level])
						dmax[level] = &(delta1[level][i]);
				for (k = 0; k < NAXS; k++)
					if (spos[k][i] > spos[k][j])
						break;            
				if (k == NAXS)
					if (delta1[level][i] > *dmin[level])
						dmin[level] = &(delta1[level][i]);
			}
        }
        store[level] = delta1[level][j];
        memcpy(lambdaold[level], lambda[level], NOBS * sizeof(double *));				
    }
	ordsavestate();
	
    do {
        order = 1;
        for (level = 1; level < NCAT; level++) {
            delta1[level][j] = gsl_ran_flat(rgen, *dmin[level], *dmax[level]);	
        }
        for (level = 1; level < NCAT; level++) {		
		    if (delta1[level-1][j] < delta1[level][j]) {
                order = 0;
            }
		}
    } while (order == 0);
	
	for (level = 1; level < NCAT; level++) {
		if (delta1[level][j] < store[level]) {
			for (i = 0; i < NOBS; i++) {
				if (lambda[level][i] == &(delta1[level][j])) {
					lambda[level][i] = &(delta0[level]);
					for (k = 0; k < steptotal; k++) {
						if (delta1[level][k] > *lambda[level][i])
							if (ordlowercorner(i, k))
								lambda[level][i] = &(delta1[level][k]);
					}
				}
			}
		}
		else if (delta1[level][j] > store[level]) {
			for (i = 0; i < NOBS; i++) {
				if (lambda[level][i] != &(delta1[level][j])) {
					if (delta1[level][j] > *lambda[level][i])
						if (ordlowercorner(i, j))
							lambda[level][i] = &(delta1[level][j]);
				}
			}
		}
	}

    mhratio = exp(ordloglik() - lold);
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        return 1;
    }
    else {
		for (level = 1; level < NCAT; level++) {
			delta1[level][j] = store[level];
            memcpy(lambda[level], lambdaold[level], NOBS * sizeof(double *));
        }
        ordrestorestate();
        scrapdvector(store);
        return 0;
    }

}

/* ************************************************************************** */

int ordupdate_delta0(int level) {

    /* Update delta0 (random walk/prior proposal Metropolis): */

    int i;
    double lold=ordloglik(), mhratio, store;
    double *dmin, *dmax;

    dmin = &(deltamin[level]);
    dmax = &(deltamax[level]);
    for (i = 0; i < steptotal; i++) {
        if (delta1[level][i] < *dmax)
            dmax = &(delta1[level][i]);
    }
    if (delta0[level-1] < *dmax)
        dmax = &(delta0[level-1]);
    if (level < NCAT-1) {
        if (delta0[level+1] > *dmin)
            dmin = &(delta0[level+1]);
    }
    ordsavestate();

    store = delta0[level];
	delta0[level] = gsl_ran_beta(rgen, 1.0 + fmin((double)steptotal, DC), 1.0) * (fmin(delta0[level] + DELTAI, *dmax) - fmax(delta0[level] - DELTAI, *dmin)) + fmax(delta0[level] - DELTAI, *dmin);

    mhratio = exp(ordloglik() - lold);
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        return 1;
    }
    else {
        delta0[level] = store;
        ordrestorestate();
        return 0;
    }
}

/* ************************************************************************** */

int ordupdate_delta0_joint() {

    /* Update delta0 jointly (random walk/prior proposal Metropolis): */

    int i, level, order;
    double lold=ordloglik(), mhratio;
	double *store;
    double **dmin, **dmax;

    store = dvector(NCAT);
    dmin = pdvector(NCAT);
    dmax = pdvector(NCAT);

    for (level = 1; level < NCAT; level++) {
        dmin[level] = &(deltamin[level]);
        dmax[level] = &(deltamax[level]);
        for (i = 0; i < steptotal; i++) {
			if (delta1[level][i] < *dmax[level])
				dmax[level] = &(delta1[level][i]);
        }
        store[level] = delta0[level];
    }
    ordsavestate();
    do {
        order = 1;
        for (level = 1; level < NCAT; level++) {
            delta0[level] = gsl_ran_flat(rgen, *dmin[level], *dmax[level]);	
        }
        for (level = 1; level < NCAT; level++) {
		    if (delta0[level-1] < delta0[level]) {
                order = 0;
            }
		}
    } while (order == 0);		

    mhratio = exp(ordloglik() - lold);

    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        return 1;
    }
    else {
		for (level = 0; level < NCAT; level++) {
			delta0[level] = store[level];
		}
        ordrestorestate();
        scrapdvector(store);		
        return 0;
    }
}

/* ************************************************************************** */

int ordupdate_beta(int j, int k, double sd) {

    /* Update beta (random walk Metropolis): */

    int i, counter=0;
    double store, lold=ordloglik(), mhratio, mean=0.0;

    store = betao[j][k];
    betao[j][k] += gsl_ran_gaussian(rgen, sd);
    memcpy(bz1old, bz1, NOBS * sizeof(double));
	
	if (GAM == 1) {
		for (i = 0; i < NOBS; i++) {
			if (include[i]) {
				mean += betao[j][(int)z[j][i]];
				counter++;
			}
		}
		mean /= counter;
		for (i = 0; i < nint[j]; i++) {
			betao[j][i] -= mean;
		}		
		for (i = 0; i < NOBS; i++) {
			if (z[j][i] == k)
				bz1[i] += (betao[j][k] - store);
			else
				bz1[i] -= mean;
		}
		if (k == 0)
			mhratio = exp(ordloglik() - lold) * gsl_ran_gaussian_pdf(betao[j][k] - (betao[j][k+1]), sqrt(sigmasq[j])) /
											 gsl_ran_gaussian_pdf(store - (betao[j][k+1]), sqrt(sigmasq[j]));
		else if (k == (nint[j] - 1))
			mhratio = exp(ordloglik() - lold) * gsl_ran_gaussian_pdf(betao[j][k] - (betao[j][k-1]), sqrt(sigmasq[j])) /
											 gsl_ran_gaussian_pdf(store - (betao[j][k-1]), sqrt(sigmasq[j]));
		else
			mhratio = exp(ordloglik() - lold) * gsl_ran_gaussian_pdf(betao[j][k] - (betao[j][k-1] + betao[j][k+1])/2.0, sqrt(sigmasq[j]/2.0)) /
											 gsl_ran_gaussian_pdf(store - (betao[j][k-1] + betao[j][k+1])/2.0, sqrt(sigmasq[j]/2.0));			
	}
	else {
		for (i = 0; i < NOBS; i++)
			bz1[i] += (betao[j][k] - store) * z[j][i];
		mhratio = exp(ordloglik() - lold);
	}	
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        return 1;
    }
    else {
		if (GAM == 1) {		
			for (i = 0; i < nint[j]; i++) {
				betao[j][i] += mean;
			}
		}
        betao[j][k] = store;
        memcpy(bz1, bz1old, NOBS * sizeof(double));
        return 0;
    }
}

/* ************************************************************************** */

void ordupdate_sigmasq(int j) {

    /* Update sigma squared (conjugate inverse gamma): */

	int k;
    double sumsq=0.0;

    for (k = 1; k < nint[j]; k++) {
		sumsq += (betao[j][k] - betao[j][k-1]) * (betao[j][k] - betao[j][k-1]);
	}
	sigmasq[j] = 1.0/gsl_ran_gamma(rgen, 0.1 + (nint[j] - 1.0)/2.0, 1.0/(0.1 + sumsq/2.0));
	return;
}

/* ************************************************************************** */

int ordupdate_alpha(int j, double sd) {

    /* Update alpha (random walk Metropolis): */

    int i;
    double store, lold=ordloglikcat(j), mhratio;

    store = alpha[j];
    alpha[j] += gsl_ran_gaussian(rgen, sd);

    memcpy(bz1old, bz1, NOBS * sizeof(double));
    for (i = 0; i < NOBS; i++)
		if (include[i])
			if (z[0][i] == j)
				bz1[i] += (alpha[j] - store);

    mhratio = exp(ordloglikcat(j) - lold) * gsl_ran_gaussian_pdf(alpha[j], sqrt(tau)) /
									     gsl_ran_gaussian_pdf(store, sqrt(tau));
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        return 1;
    }
    else {
        alpha[j] = store;
        memcpy(bz1, bz1old, NOBS * sizeof(double));
        return 0;
    }
}

/* ************************************************************************** */

void ordupdate_tau() {

    /* Update tau (conjugate inverse gamma): */

	int j;
    double sumsq=0.0;

    for (j = 0; j < NCLUSTER; j++) {
		sumsq += alpha[j] * alpha[j];
	}
	tau = 1.0/gsl_ran_gamma(rgen, 0.1 + NCLUSTER/2.0, 1.0/(0.1 + sumsq/2.0));
	return;
}

/* ************************************************************************** */

void ordupdate_dimtotals() {
	int i, j;

    for (i = 0; i < NAXS; i++) {
		dimtotals[i] = 0;
		for (j = 0; j < NPPS; j++) {
			if (settozero[j][i] == 0)
				dimtotals[i] += steps[j];
		}
	}
	return;
}

/* ************************************************************************** */

int ordbirth(int pid) {
    /* Propose new changepoint: */
    int i, j, level, order;
    double lold=ordloglik(), mhratio;
    double **dmin, **dmax;

    dmin = pdvector(NCAT);
    dmax = pdvector(NCAT);

    pp[steptotal] = pid;

    for (i = 0; i < NAXS; i++) {
		if ((dimtotals[i] == 0) && (settozero[pid][i] == 0)) {
			if (gsl_rng_uniform(rgen) < INVPROB) {
				direction[i] = 1;
			}
			else {
				direction[i] = -1;				
				ordinvert(i);
			}
		}
        if (settozero[pid][i])
            spos[i][steptotal] = 0.0;
        else {
            spos[i][steptotal] = gsl_ran_flat(rgen, smin[i], smax[i]);
        }
    }
    ordsavestate();
    delta1[0][steptotal] = DUPPER;
    for (level = 1; level < NCAT; level++) {
        dmin[level] = &(delta0[level]);
        dmax[level] = &(deltamax[level]);
        for (i = 0; i < steptotal; i++) {
            for (j = 0; j < NAXS; j++)
                if (spos[j][i] < spos[j][steptotal])
                    break;            
            if (j == NAXS) 
                if (delta1[level][i] < *dmax[level])
                    dmax[level] = &(delta1[level][i]);
            for (j = 0; j < NAXS; j++)
                if (spos[j][i] > spos[j][steptotal])
                    break;            
            if (j == NAXS)
                if (delta1[level][i] > *dmin[level])
                    dmin[level] = &(delta1[level][i]);
        }
        memcpy(lambdaold[level], lambda[level], NOBS * sizeof(double *));
    }
    do {
        order = 1;
        for (level = 1; level < NCAT; level++) {
            delta1[level][steptotal] = gsl_ran_flat(rgen, *dmin[level], *dmax[level]);	
        }
        for (level = 1; level < NCAT; level++) {		
		    if (delta1[level-1][steptotal] < delta1[level][steptotal]) {
                order = 0;
            }
		}
    } while (order == 0);
    for (level = 0; level < NCAT; level++) {
        for (i = 0; i < NOBS; i++) {
			if (delta1[level][steptotal] > *lambda[level][i])
				if (ordlowercorner(i, steptotal))
					lambda[level][i] = &(delta1[level][steptotal]);
        }
    }
    steptotal++;
    steps[pid]++;
	ordupdate_dimtotals();

    scrappdvector(dmin);
    scrappdvector(dmax);

    /* Accept/reject: */
    mhratio = exp(ordloglik() - lold) * rho[pid] * vol[pid] / steps[pid];
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        return 1;
    }
    else {
        for (level = 1; level < NCAT; level++) {
            memcpy(lambda[level], lambdaold[level], NOBS * sizeof(double *));
        }
        steptotal--;
        steps[pid]--;
		ordupdate_dimtotals();
		for (i = 0; i < NAXS; i++) {
			if ((dimtotals[i] == 0) && (settozero[pid][i] == 0)) {
				if (direction[i] == -1)
					ordinvert(i);
				direction[i] = 0;
			}
		}
        ordrestorestate();
        return 0;
    }
}

/* ************************************************************************** */

void ordkill(int j) {
    /* Kill an existing changepoint: */
    int i=0, k, l, level;
    for (level = 0; level < NCAT; level++) {
        for (i = 0; i < NOBS; i++) {
			if (lambda[level][i] == &(delta1[level][j])) {
				lambda[level][i] = &(delta0[level]);
				for (k = 0; k < steptotal; k++) {
					if (k != j) {
						if (delta1[level][k] > *lambda[level][i])
							if (ordlowercorner(i, k))
								lambda[level][i] = &(delta1[level][k]);
					}
				}
			}
			if (lambda[level][i] != &(delta0[level])) {
				if (lambda[level][i] > &(delta1[level][j]))
					lambda[level][i]--;
			}
        }
    }
    steptotal--;
    steps[pp[j]]--;
	ordupdate_dimtotals();	
    for (i = j; i < steptotal; i++) {
        for (l = 0; l < NAXS; l++)
            spos[l][i] = spos[l][i+1];
        for (level = 0; level < NCAT; level++)
            delta1[level][i] = delta1[level][i+1];
        pp[i] = pp[i+1];
    }
	return;
}

/* ************************************************************************** */

void ordgc() {
	int i, j, sum=0, total=steptotal;
	int *count;
    count = ivector(MAXTOTAL);
    for (j = 0; j < total; j++) {
		count[j] = 0;
	}	
	for (i = 0; i < NOBS; i++) {
		if (include[i]) {
			if (d[i] == NCAT-1) {
				if (lambda[d[i]][i] != &(delta0[d[i]]))
					count[lambda[d[i]][i] - &(delta1[d[i]][0])]++;
			}
			else if (d[i] == 0) {
				if (lambda[d[i]+1][i] != &(delta0[d[i]+1]))			
					count[lambda[d[i]+1][i] - &(delta1[d[i]+1][0])]++;
			}
			else {
				if (lambda[d[i]][i] != &(delta0[d[i]])) {					
					count[lambda[d[i]][i] - &(delta1[d[i]][0])]++;
				}
				if (lambda[d[i]+1][i] != &(delta0[d[i]+1])) {					
					count[lambda[d[i]+1][i] - &(delta1[d[i]+1][0])]++;
				}
			}
		}
	}
    for (j = 0; j < total; j++) {
		if (count[j] == 0) { 
			for (i = 0; i < NAXS; i++)
				if (spos[i][j-sum] > 0.25)
					break;
			if (i == NAXS) {
				ordkill(j-sum);
				sum++;
			}
		}
	}
    scrapivector(count);	
	return;
}

/* ************************************************************************** */

int orddeath(int pid) {
    /* Propose to remove an existing changepoint: */
    int i=0, j=gsl_rng_uniform_int(rgen, steps[pid]), k=-1, l, level;
	int *olddirection;
    double lold=ordloglik(), mhratio;
    double *sold, *store;
	olddirection = ivector(NAXS);	
    sold = dvector(NAXS);
    store = dvector(NCAT);

    while (i <= j) {
        k++;
        if (pp[k] == pid)
            i++;
    }
    j = k;

    ordsavestate();
    for (i = 0; i < NAXS; i++) {
        sold[i] = spos[i][j];
    }
    for (level = 0; level < NCAT; level++) {
        store[level] = delta1[level][j];
        memcpy(lambdaold[level], lambda[level], NOBS * sizeof(double *));
        for (i = 0; i < NOBS; i++) {
			if (lambda[level][i] == &(delta1[level][j])) {
				lambda[level][i] = &(delta0[level]);
				for (k = 0; k < steptotal; k++) {
					if (k != j) {
						if (delta1[level][k] > *lambda[level][i])
							if (ordlowercorner(i, k))
								lambda[level][i] = &(delta1[level][k]);
					}
				}
			}
			if (lambda[level][i] != &(delta0[level])) {
				if (lambda[level][i] > &(delta1[level][j]))
					lambda[level][i]--;
			}
        }
    }
    steptotal--;
    steps[pid]--;
	ordupdate_dimtotals();
	for (i = 0; i < NAXS; i++) {
		olddirection[i] = direction[i];
		if ((dimtotals[i] == 0) && (settozero[pid][i] == 0)) {
			if (direction[i] == -1)
				ordinvert(i);
			direction[i] = 0;
		}
	}	
    for (i = j; i < steptotal; i++) {
        for (l = 0; l < NAXS; l++)
            spos[l][i] = spos[l][i+1];
        for (level = 0; level < NCAT; level++)
            delta1[level][i] = delta1[level][i+1];
        pp[i] = pp[i+1];
    }

    /* Accept/reject: */
    mhratio = exp(ordloglik() - lold) / rho[pid] * (steps[pid] + 1) / vol[pid];
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        scrapdvector(sold);
        scrapdvector(store);
        return 1;
    }
    else {
        for (level = 0; level < NCAT; level++) {
            memcpy(lambda[level], lambdaold[level], NOBS * sizeof(double *));
        }
        for (i = steptotal; i > j; i--) {
            for (l = 0; l < NAXS; l++)
                spos[l][i] = spos[l][i-1];
            for (level = 0; level < NCAT; level++)
                delta1[level][i] = delta1[level][i-1];
            pp[i] = pp[i-1];
        }
        for (l = 0; l < NAXS; l++)
            spos[l][j] = sold[l];
        for (level = 0; level < NCAT; level++)
            delta1[level][j] = store[level];
        pp[j] = pid;
        steptotal++;
        steps[pid]++;
		ordupdate_dimtotals();		
		for (i = 0; i < NAXS; i++) {
			if ((olddirection[i] == -1) && (direction[i] == 0))
				ordinvert(i);
			direction[i] = olddirection[i];
		}		
        ordrestorestate();
        scrapdvector(sold);
        scrapdvector(store);
        return 0;
    }
}

/* ************************************************************************** */

int ordmove(int pid) {
    /* Propose to move an existing changepoint: */
    int i=0, j = gsl_rng_uniform_int(rgen, steps[pid]), k=-1, l, level;

    while (i <= j) {
        k++;
        if (pp[k] == pid)
            i++;
    }
    j = k;

    double lold=ordloglik(), mhratio;
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
        if (i != j) {
            for (l = 0; l < NAXS; l++) {
                if (spos[l][i] <= spos[l][j] && spos[l][i] > *lowerlim[l])
                    lowerlim[l] = &(spos[l][i]);
                if (spos[l][i] >= spos[l][j] && spos[l][i] < *upperlim[l])
                    upperlim[l] = &(spos[l][i]);
            }
        }
    }
    for (i = 0; i < NAXS; i++) {
        if (settozero[pid][i])
            spos[i][j] = 0.0;
        else
            spos[i][j] = gsl_ran_flat(rgen, *lowerlim[i], *upperlim[i]);
    }
    ordsavestate();
    for (level = 0; level < NCAT; level++) {
        memcpy(lambdaold[level], lambda[level], NOBS * sizeof(double *));
        for (i = 0; i < NOBS; i++) {
			if (lambda[level][i] == &(delta1[level][j])) {
				lambda[level][i] = &(delta0[level]);
				for (k = 0; k < steptotal; k++) {
						if (delta1[level][k] > *lambda[level][i])
							if (ordlowercorner(i, k))
								lambda[level][i] = &(delta1[level][k]);
				}
			}
			else if (lambda[level][i] != &(delta1[level][j])) {
				if (delta1[level][j] > *lambda[level][i])
					if (ordlowercorner(i, j))
						lambda[level][i] = &(delta1[level][j]);
			}
        }
    }

    /* Accept/reject: */
    mhratio = exp(ordloglik() - lold);
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        scrapdvector(sold);
        scrappdvector(lowerlim);
        scrappdvector(upperlim);
        return 1;
    }
    else {
        for (level = 0; level < NCAT; level++)
            memcpy(lambda[level], lambdaold[level], NOBS * sizeof(double *));
        for (i = 0; i < NAXS; i++)
            spos[i][j] = sold[i];
        ordrestorestate();
        scrapdvector(sold);
        scrappdvector(lowerlim);
        scrappdvector(upperlim);
        return 0;
    }
}

/* ************************************************************************** */

int orddeath_birth(int pid) {
    /* Combined death-birth move: */
    int i=0, j=gsl_rng_uniform_int(rgen, steps[pid]), k=-1, newpid, l, level, order;
	int *olddirection;
    double lold=ordloglik(), mhratio;
    double *sold, *store;
    double **dmin, **dmax;
	olddirection = ivector(NAXS);
    sold = dvector(NAXS);
    store = dvector(NCAT);
    dmin = pdvector(NCAT);
    dmax = pdvector(NCAT);

    while (i <= j) {
        k++;
        if (pp[k] == pid)
            i++;
    }
    j = k;
    newpid = pps[gsl_rng_uniform_int(rgen, NPPS)];
    /* newpid = gsl_ran_discrete(rgen, discd); */

    /* Store old values: */

    ordsavestate();
    for (level = 0; level < NCAT; level++) {
        store[level] = delta1[level][j];
        memcpy(lambdaold[level], lambda[level], NOBS * sizeof(double *));
    }
    for (i = 0; i < NAXS; i++) {
        sold[i] = spos[i][j];
    }

    /* Remove pointers to the removed point: */

    for (level = 0; level < NCAT; level++) {
        for (i = 0; i < NOBS; i++) {
			if (lambda[level][i] == &(delta1[level][j])) {
				lambda[level][i] = &(delta0[level]);
				for (k = 0; k < steptotal; k++) {
					if (k != j) {
						if (delta1[level][k] > *lambda[level][i])
							if (ordlowercorner(i, k))
								lambda[level][i] = &(delta1[level][k]);
					}
				}
			}
        }
    }
    steps[pid]--;
	ordupdate_dimtotals();
	for (i = 0; i < NAXS; i++) {
		olddirection[i] = direction[i];
		if ((dimtotals[i] == 0) && (settozero[pid][i] == 0)) {
			if (direction[i] == -1)
				ordinvert(i);
			direction[i] = 0;
		}
	}	

    /* New point: */

    pp[j] = newpid;
    for (i = 0; i < NAXS; i++) {
		if ((dimtotals[i] == 0) && (settozero[newpid][i] == 0)) {
			if (gsl_rng_uniform(rgen) < INVPROB) {
				direction[i] = 1;
			}
			else {
				direction[i] = -1;				
				ordinvert(i);
			}
		}		
        if (settozero[newpid][i])
            spos[i][j] = 0.0;
        else {
            spos[i][j] = gsl_ran_flat(rgen, smin[i], smax[i]);
        }
    }

    delta1[0][j] = DUPPER;
    for (level = 1; level < NCAT; level++) {
        dmin[level] = &(delta0[level]);
        dmax[level] = &(deltamax[level]);
        for (i = 0; i < steptotal; i++) {
            if (i != j) {
                for (l = 0; l < NAXS; l++)
                    if (spos[l][i] < spos[l][j])
                        break;
                if (l == NAXS)
                    if (delta1[level][i] < *dmax[level])
                        dmax[level] = &(delta1[level][i]);
                for (l = 0; l < NAXS; l++)
                    if (spos[l][i] > spos[l][j])
                        break;
                if (l == NAXS)
                    if (delta1[level][i] > *dmin[level])
                        dmin[level] = &(delta1[level][i]);
            }
        }
    }
    do {
        order = 1;
        for (level = 1; level < NCAT; level++) {
            delta1[level][j] = gsl_ran_flat(rgen, *dmin[level], *dmax[level]);
            if (delta1[level-1][j] < delta1[level][j]) {
                order = 0;
                break;
            }
        }
    } while (order == 0);
    for (level = 0; level < NCAT; level++) {
        for (i = 0; i < NOBS; i++) {
			if (delta1[level][j] > *lambda[level][i])
				if (ordlowercorner(i, j))
					lambda[level][i] = &(delta1[level][j]);
        }
    }
    steps[newpid]++;
	ordupdate_dimtotals();

    /* Accept/reject: */
    mhratio = exp(ordloglik() - lold) * (rho[newpid] / rho[pid]) *
              ((double)(steps[pid] + 1) / steps[newpid]) * (vol[newpid] / vol[pid]);
    if (gsl_rng_uniform_pos(rgen) < GSL_MIN_DBL(1.0, mhratio)) {
        scrapdvector(sold);
        scrapdvector(store);
        scrappdvector(dmin);
        scrappdvector(dmax);
        return 1;
    }
    else {
        for (level = 0; level < NCAT; level++)
            memcpy(lambda[level], lambdaold[level], NOBS * sizeof(double *));
        steps[newpid]--;
        steps[pid]++;
        pp[j] = pid;
        for (i = 0; i < NAXS; i++)
            spos[i][j] = sold[i];
        for (level = 0; level < NCAT; level++)
            delta1[level][j] = store[level];
		ordupdate_dimtotals();
		
		for (i = 0; i < NAXS; i++) {
			if ((olddirection[i] == 0) && (direction[i] == -1))
				ordinvert(i);
			if ((olddirection[i] == -1) && (direction[i] == 1))
				ordinvert(i);
			if ((olddirection[i] == 1) && (direction[i] == -1))
				ordinvert(i);
			if ((olddirection[i] == -1) && (direction[i] == 0))
				ordinvert(i);			
			direction[i] = olddirection[i];
		}
        ordrestorestate();
        scrapivector(olddirection);
        scrapdvector(sold);
        scrapdvector(store);
        scrappdvector(dmin);
        scrappdvector(dmax);
        return 0;
    }
}

/* ************************************************************************** */

int ordcheckpartialordering() {
    int i, j, k, level;
    for (level = 0; level < NCAT; level++) {
        for (i = 0; i < steptotal; i++) {
            if ((delta1[level][i] < deltamin[level]) || (delta1[level][i] > deltamax[level])) {
                Rprintf("Partial ordering violated (deltamin/deltamax).\n");
                return(1);
            }
            if (delta0[level] > delta1[level][i]) {
                Rprintf("Partial ordering violated (delta0).\n");
                return(1);
            }
            for (j = 0; j < steptotal; j++) {
                if (j != i) {
                    for (k = 0; k < NAXS; k++)
                        if (spos[k][j] > spos[k][i])
                            break;
                    if (k == NAXS)
                        if (delta1[level][j] > delta1[level][i])
                            goto error1;
                    for (k = 0; k < NAXS; k++)
                        if (spos[k][j] < spos[k][i])
                            break;
                    if (k == NAXS)
                        if (delta1[level][j] < delta1[level][i])
                            goto error1;
                }
            }
            if (level > 0)
                if (delta1[level-1][i] < delta1[level][i])
                    goto error2;
            if (level < NCAT-1)
                if (delta1[level+1][i] > delta1[level][i])
                    goto error2;
        }
    }
    return 0;
    error1:
    Rprintf("Partial ordering violated (levels).\n");
    error2:
    Rprintf("Partial ordering violated (survival).\n");
    return 1;
}

/* ************************************************************************** */

void ordsampler(int *iargs, double *dargs, int *idata, double *ddata,
                int *isettozero, int *steptotalf, int *stepsf,
                double *rhof, double *likf, double *tauf, double *alphaf, double *betaf, double *lambdaf, double *sigmasqf) {

    /* Copy constant arguments: */

    NITER = iargs[0];
    BURNIN = iargs[1];
    ADAPT = iargs[2];
    REFRESH = iargs[3];
    THIN = iargs[4];
    NOBS = iargs[5];
    NCOV = iargs[6];
    NAXS = iargs[7];
    BIRTHDEATH = iargs[8];
    NPPS = iargs[9];
    NCAT = iargs[10];
    NCLUSTER = iargs[11];
    LOGIT = iargs[12];
    SEED = iargs[13];
	GAM = iargs[14];

    RHOA = dargs[0];
    RHOB = dargs[1];
    DELTAI = dargs[2];
    DLOWER = dargs[3];
    DUPPER = dargs[4];
	INVPROB = dargs[5];
	DC = dargs[6];

    /* Local variables: */

    int i, j, k, l, m, nupdates, counter, scounter=0, nsave=(NITER - BURNIN)/THIN, delta0jointprop, delta0jointacc, deltajointprop, deltajointacc, nmax;
    long int *deltaprop, *deltaacc, *delta0prop, *delta0acc;
    long int **prop, **acc;
    double *size;
    double sum;

    deltaprop = livector(NCAT);
    deltaacc = livector(NCAT);
    delta0prop = livector(NCAT);
    delta0acc = livector(NCAT);
    prop = limatrix(NPPS, 4);
    acc = limatrix(NPPS, 4);
    size = dvector(NPPS);

    gsl_set_error_handler_off();

    /* Counters to zero: */

	delta0jointprop=0;
	delta0jointacc=0;
	deltajointprop=0;
	deltajointacc=0;
    for (k = 0; k < NCAT; k++) {
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
    direction = ivector(NAXS);
    dimtotals = ivector(NAXS);
	nint = ivector(NCOV);
    
    settozero = imatrix(NPPS, NAXS);

    alphacnt = limatrix(1, 2);

    smin = dvector(NAXS);
    smax = dvector(NAXS);
    vol = dvector(NPPS);
    delta0 = dvector(NCAT);
    delta0old = dvector(NCAT);
    deltamin = dvector(NCAT);
    deltamax = dvector(NCAT);
    deltaminold = dvector(NCAT);
    deltamaxold = dvector(NCAT);
    rho = dvector(NPPS);
    alpha = dvector(NCLUSTER);
    bz1 = dvector(NOBS);
    bz1old = dvector(NOBS);
    sigmasq = dvector(NCOV);

    z = dmatrix(NCOV, NOBS);
    x = dmatrix(NAXS, NOBS);
    spos = dmatrix(NAXS, MAXTOTAL);
    alphasd = dmatrix(1, 3);
    delta1 = dmatrix(NCAT, MAXTOTAL);
    delta1old = dmatrix(NCAT, MAXTOTAL);
	
    lambda = pdmatrix(NCAT, NOBS);
    lambdaold = pdmatrix(NCAT, NOBS);

    rgen = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rgen, SEED);

    /* Read data: */

    Rprintf("Read data:\n");
    for (j = 0; j < NCOV; j++) {
		nint[j] = 1;		
	}
    for (i = 0; i < NOBS; i++) {
        predict[i] = idata[vidx(i, 0, NOBS)];
        include[i] = idata[vidx(i, 1, NOBS)];
        d[i] = idata[vidx(i, 2, NOBS)];
        if (i < 5)
            Rprintf("%d %d %d ", predict[i], include[i], d[i]);

        counter = 0;
        for (j = 0; j < NAXS; j++) {
            x[j][i] = ddata[vidx(i, counter, NOBS)];
            counter++;
            if (i < 5)
                Rprintf("%f ", x[j][i]);
        }
        for (j = 0; j < NCOV; j++) {
            z[j][i] = ddata[vidx(i, counter, NOBS)];
			if (GAM == 1 && j > 0)
				nint[j] = (((int)z[j][i]+1) > nint[j]) ? ((int)z[j][i]+1) : nint[j];
            counter++;
            if (i < 5)
                Rprintf("%f ", z[j][i]);
        }
        if (i < 5)
            Rprintf("\n");
    }
	nmax = 0;
    for (j = 0; j < NCOV; j++) {
		nmax = (nint[j] > nmax) ? nint[j] : nmax;
		if (GAM==1) {
			Rprintf("Number of categories for x[%d]: %d", j, nint[j]);
			Rprintf("\n");
		}
	}
	if (GAM==1) {
		Rprintf("Max number of categories %d", nmax);
		Rprintf("\n");
	}
    betacnt = limatrix(NCOV * nmax, 2);
    betasd = dmatrix(NCOV * nmax, 3);
    betao = dmatrix(NCOV, nmax);

    Rprintf("Point process configuration:\n");
    for (i = 0; i < NPPS; i++) {
        pps[i] = i;
        for (j = 0; j < NAXS; j++) {
            settozero[i][j] = isettozero[vidx(i, j, NPPS)];
            Rprintf("%d ", settozero[i][j]);
        }
        Rprintf("\n");
    }
    Rprintf("Number of point processes: %d\n", NPPS);
    Rprintf("Number of covariate axes: %d\n", NAXS);
    Rprintf("Number of outcome categories: %d\n", NCAT);
    Rprintf("Number of observations: %d\n", NOBS);

    /* Set initial values: */

    for (i = 0; i < NAXS; i++) {
		direction[i] = 0;
		dimtotals[i] = 0;
        smin[i] = 0.0;
        smax[i] = 1.0;
    }
    steptotal = 0;
    for (i = 0; i < NPPS; i++) {
        vol[i] = 1.0;
        steps[i] = 0;
        rho[i] = RHOA/RHOB;
        maxstep[i] = NOBS;
    }
	for (i = 0; i < NOBS; i++) {
		if (include[i])
			ninc++;
	}	
    for (j = 0; j < NCAT; j++) {
		counter = 0;
		for (i = 0; i < NOBS; i++) {
			if (include[i])
				counter += (d[i] >= j);
		}
		if (LOGIT)
			delta0[j] = (j == 0) ? DUPPER : log(((double)counter/ninc)/(1.0 - (double)counter/ninc));
		else 
			delta0[j] = (double)counter/ninc;			
        Rprintf("Init. delta0[%d] = %f\n", j, delta0[j]);		
        deltamin[j] = DLOWER;
        deltamax[j] = DUPPER;
    }

    delta0[0] = DUPPER;
	tau=0.1;
    for (i = 0; i < NCLUSTER; i++) {
        alpha[i] = 0.0;
	}
	setupproposal(alphacnt[0], alphasd[0], 0.1);
	
    for (i = 0; i < NCOV; i++) {
		for (j = 0; j < nint[i]; j++) {
			betao[i][j] = 0.0;
			m = vidx(i, j, NCOV);
			setupproposal(betacnt[m], betasd[m], 0.1);
		}
		sigmasq[i] = 0.1;
    }

    for (i = 0; i < NOBS; i++) {
        bz1[i] = 0.0;
    }
    for (i = 0; i < NOBS; i++) {
        for (j = 0; j < NCAT; j++) {
            lambda[j][i] = &(delta0[j]);
        }
    }

    /* Gibbs sampler loop: */

    time_t t0, t1, t2;
    t0 = time(NULL);
    t1 = time(NULL);
    double tdiff;

    Rprintf("Starting Gibbs sampler:\n");

    for (i = 1; i <= NITER; i++) {

        /* Rprintf("Adjust proposals:\n"); */

        /* Adjust proposals: */

        if ((i > 100 && i < ADAPT && i % 100 == 0) || i == ADAPT) {
            adjustproposal(alphacnt[0], alphasd[0], i, ADAPT, 0.1);
            for (k = 1; k < NCOV; k++) {
				for (j = 0; j < nint[k]; j++) {
					m = vidx(k, j, NCOV);
					adjustproposal(betacnt[m], betasd[m], i, ADAPT, 0.1);
				}
			}
        }

        /* Rprintf("Birth/death moves:\n"); */

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
                    acc[l][0] += ordbirth(l);
                }
            }
            else if (steps[l] == maxstep[l]) {
                if (gsl_rng_uniform(rgen) < 0.5) {
                    prop[l][1]++;
                    acc[l][1] += orddeath(l);
                }
            }
            else {
                if (gsl_rng_uniform(rgen) < 0.5) {
                    prop[l][0]++;
                    acc[l][0] += ordbirth(l);
                }
                else {
                    prop[l][1]++;
                    acc[l][1] += orddeath(l);
                }
            }
        }

        /* Rprintf("Combined death/birth moves:\n"); */

        /* Combined death/birth moves: */

        for (k = 0; k < BIRTHDEATH; k++) {
            l = pps[gsl_rng_uniform_int(rgen, NPPS)];
            l = gsl_ran_discrete(rgen, discd);
            if (steps[l] > 0) {
                prop[l][2]++;
                acc[l][2] += orddeath_birth(l);
            }
        }

        /* Rprintf("Position change moves:\n"); */

        /* Position change moves: */

        nupdates = GSL_MIN_INT(steptotal, 100);
        for (k = 0; k < nupdates; k++) {
            l = pp[gsl_rng_uniform_int(rgen, steptotal)];
            prop[l][3]++;
            acc[l][3] += ordmove(l);
        }

		/* ordgc(); */

        /* Rprintf("Rho updates:\n"); */

        for (k = 0; k < NPPS; k++) {
            /* l = pps[gsl_rng_uniform_int(rgen, NPPS)]; */
			/* if (l == NPPS - 1) */
			ordupdate_rho(k);
        }

        /* Other parameters: */

        /* Rprintf("Delta0 updates:\n"); */

		if (!(DC > 0.0)) {
			delta0jointacc += ordupdate_delta0_joint();
			delta0jointprop++;
		}
        for (k = 0; k < NCAT; k++) {
            l = gsl_rng_uniform_int(rgen, NCAT-1) + 1;
            delta0acc[l] += ordupdate_delta0(l);
            delta0prop[l]++;
        }

        /* Rprintf("Delta updates:\n"); */

        nupdates = GSL_MIN_INT(steptotal, 100);
        for (k = 0; k < nupdates; k++) {
            deltajointacc += ordupdate_delta_joint(gsl_rng_uniform_int(rgen, steptotal));
            deltajointprop++;
        }

        nupdates = GSL_MIN_INT(steptotal * (NCAT-1), 100);
        for (k = 0; k < nupdates; k++) {
            l = gsl_rng_uniform_int(rgen, NCAT-1) + 1;
            deltaacc[l] += ordupdate_delta(gsl_rng_uniform_int(rgen, steptotal), l);
            deltaprop[l]++;
        }

        /* Rprintf("Alpha updates:\n"); */

		if (i > round(ADAPT/5)) {			
			if (NCLUSTER > 1) {
				for (k = 0; k < NCLUSTER; k++) {
					alphacnt[0][0] += ordupdate_alpha(k, alphasd[0][0]);
					alphacnt[0][1]++;
				}
				ordupdate_tau();
			}

			/* Rprintf("Beta updates:\n"); */

			for (j = 1; j < NCOV; j++) {
				if (GAM == 1) {
					for (k = 0; k < nint[j]; k++) {
						l = gsl_rng_uniform_int(rgen, nint[j]);
						m = vidx(j, l, NCOV);
						betacnt[m][0] += ordupdate_beta(j, l, betasd[m][0]);
						betacnt[m][1]++;
					}
					ordupdate_sigmasq(j);
				}
				else {
					betacnt[j][0] += ordupdate_beta(j, 0, betasd[j][0]);
					betacnt[j][1]++;
				}
			}
		}
		
        /* Rprintf("Write the output:\n"); */

        /* Write the output: */

        if (i > BURNIN && i % THIN == 0) {
			if (1) {
				/* Save predictions: */
				for (k = 0; k < NCAT; k++) {
					lambdaf[vidx(scounter * NCAT + k, 0, nsave * NCAT)] = k;
					counter = 0;
					if (LOGIT) {
						for (j = 0; j < NOBS; j++) {
							if (predict[j]) {
								lambdaf[vidx(scounter * NCAT + k, counter + 1, nsave * NCAT)] = (k == 0) ? 1.0 : 1.0/(1.0 + exp(-(*lambda[k][j] + bz1[j])));
								counter++;
							}
						}
					}
					else {
						for (j = 0; j < NOBS; j++) {
							if (predict[j]) {
								lambdaf[vidx(scounter * NCAT + k, counter + 1, nsave * NCAT)] = *lambda[k][j];
								counter++;
							}
						}					
					}
				}
			}
			if (0) {
				/* Calculate running averages: */
				for (k = 0; k < NCAT; k++) {
					lambdaf[vidx(0 * NCAT + k, 0, 1 * NCAT)] = k;
					counter = 0;
					if (LOGIT) {
						for (j = 0; j < NOBS; j++) {
							if (predict[j]) {
								lambdaf[vidx(0 * NCAT + k, counter + 1, 1 * NCAT)] = (k == 0) ? 1.0 : (lambdaf[vidx(0 * NCAT + k, counter + 1, 1 * NCAT)] * scounter + 1.0/(1.0 + exp(-(*lambda[k][j] + bz1[j]))))/(scounter+1);
								counter++;
							}
						}
					}
					else {
						for (j = 0; j < NOBS; j++) {
							if (predict[j]) {
								lambdaf[vidx(0 * NCAT + k, counter + 1, 1 * NCAT)] = (lambdaf[vidx(0 * NCAT + k, counter + 1, 1 * NCAT)] * scounter + *lambda[k][j])/(scounter+1);
								counter++;
							}
						}					
					}
				}
			}
            for (k = 0; k < NPPS; k++) {
                rhof[vidx(scounter, k, nsave)] = rho[k];
                stepsf[vidx(scounter, k, nsave)] = steps[k];
            }
            for (k = 0; k < NCLUSTER; k++)
                alphaf[vidx(scounter, k, nsave)] = alpha[k];			
            for (k = 0; k < NCOV; k++) {
				for (l = 0; l < nint[k]; l++)
					betaf[vidx(scounter, vidx(k, l, NCOV), nsave)] = betao[k][l];
				sigmasqf[vidx(scounter, k, nsave)] = sigmasq[k];				
			}

            steptotalf[scounter] = steptotal;
            likf[scounter] = ordloglik();
            tauf[scounter] = tau;
            scounter++;   			
        }

        /* printf("Checking partial ordering: "); */
        if (ordcheckpartialordering())
            goto closeoutputfiles;
        /* printf("Done.\n"); */

        /* Time interval: */

        if (i % REFRESH == 0) {
			Rprintf("Acceptance ratio for delta0 joint updates: %f\n", (double)delta0jointacc/delta0jointprop);
			Rprintf("Acceptance ratio for delta joint updates: %f\n", (double)deltajointacc/deltajointprop);
            for (k = 0; k < NCAT; k++) {
                Rprintf("Acceptance ratio for delta0[%d] (=%f) updates: %f\n", k, delta0[k], (double)delta0acc[k]/delta0prop[k]);
            }
            for (l = 0; l < NCAT; l++) {
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
			for (l = 0; l < NPPS; l++)
                Rprintf("n[%d]=%d, rho[%d]=%f\n", l, steps[l], l, rho[l]);
            t2 = time(NULL);
            tdiff = difftime(t2, t1);
            Rprintf("%d iterations (%d s)\n", i, (int)tdiff);
            Rprintf("%d points, l=%f\n", steptotal, ordloglik());
            Rprintf("%f ", tau);			
            Rprintf("\n");
            for (k = 0; k < GSL_MIN(NCLUSTER, 5); k++)
                Rprintf("%f ", alpha[k]);
            Rprintf("\n");
			
			for (k = 0; k < NCOV; k++) {
				for (l = 0; l < nint[k]; l++)
					Rprintf("%f ", betao[k][l]);
				Rprintf("\n");				
			}
            Rprintf("\n");
            for (k = 0; k < NCOV; k++) {
				Rprintf("%f ", sigmasq[k]);			
			}
            Rprintf("\n");
            for (k = 0; k < NAXS; k++) {
				Rprintf(" %d: %d (%d)", k+1, direction[k], dimtotals[k]);
			}
            Rprintf("\n");
			
			/*
			for (l = 0; l < 5; l++) {	
				for (k = 0; k < NAXS; k++) {
						Rprintf("%f ", x[k][l]);
				}
				Rprintf("\n");
			}
			*/
			
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
    Rprintf("Acceptance ratio for delta0 joint updates: %f\n", (double)delta0jointacc/delta0jointprop);
    Rprintf("Acceptance ratio for delta joint updates: %f\n", (double)deltajointacc/deltajointprop);
    for (k = 0; k < NCAT; k++) {
        Rprintf("Acceptance ratio for delta0[%d] updates: %f\n", k, (double)delta0acc[k]/delta0prop[k]);
        Rprintf("Acceptance ratio for delta[%d] updates: %f\n", k, (double)deltaacc[k]/deltaprop[k]);
    }
    Rprintf("Acceptance ratio for alpha updates (sd): %f (%f)\n", (double)alphacnt[0][0]/alphacnt[0][1], alphasd[0][0]);
    for (k = 0; k < NCOV; k++) {
		for (l = 0; l < nint[k]; l++) {
			m = vidx(k, l, NCOV);
			Rprintf("Acceptance ratio for beta[%d][%d] updates (sd): %f (%f)\n", k, l, (double)betacnt[m][0]/betacnt[m][1], betasd[m][0]);
		}
	}

    /* Free memory: */

    closeoutputfiles:

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
    scrapivector(direction);
    scrapivector(dimtotals);
	scrapivector(nint);    
	
    scrapimatrix(settozero, NPPS);
    scraplimatrix(betacnt, NCOV * nmax);

    scrapdvector(smin);
    scrapdvector(smax);
    scrapdvector(vol);
    scrapdvector(delta0);
    scrapdvector(delta0old);
    scrapdvector(deltamin);
    scrapdvector(deltamax);
    scrapdvector(deltaminold);
    scrapdvector(deltamaxold);
    scrapdvector(rho);
    scrapdvector(alpha);
    scrapdvector(bz1);
    scrapdvector(bz1old);
    scrapdvector(sigmasq);

    scrapdmatrix(z, NCOV);
    scrapdmatrix(x, NAXS);
    scrapdmatrix(spos, NAXS);
    scrapdmatrix(betasd, NCOV * nmax);
    scrapdmatrix(delta1, NCAT);
    scrapdmatrix(delta1old, NCAT);
    scrapdmatrix(betao, NCOV);

    scrappdmatrix(lambda, NCAT);
    scrappdmatrix(lambdaold, NCAT);

    gsl_rng_free(rgen);
    return;
}

/* ************************************************************************** */
