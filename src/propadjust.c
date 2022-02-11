/* ************************************************************************** */
/* File:        propadjust.c                                                  */
/* Author:      Olli Saarela (olli.saarela@utoronto.ca)                       */
/* Summary:     Functions for adjusting Metropolis proposal etc.              */
/* ************************************************************************** */

#include <stdio.h>
#include <math.h>

#include "propadjust.h"

/* ************************************************************************** */

double sdadj(double sd, double rate) {

    /* Adjust proposal standard deviation: */

    if (rate < 0.5)
        sd *= 0.618 + 0.382*(rate/0.5);
    else if (rate > 0.5)
        sd *= 1.0 + 0.618*((rate-0.5)/0.5);

    return sd;
}

/* ************************************************************************** */

void adjustproposal(long int *counter, double *sd, int iter, int adapt, double init) {

    /* Adjust a Metropolis updater: */

    double rate, weight;

    if (counter[1] > 0 && iter <= adapt) {
        rate = (double)counter[0]/counter[1];
        weight = 1.0 - 2.0*(fabs(rate - 0.5));
        if ((iter >= round(adapt/5)) & (iter <= adapt)) {
            sd[1] += weight * sd[0];
            sd[2] += weight;
        }
        sd[0] = sdadj(sd[0], rate);
        counter[0] = 0;
        counter[1] = 0;
    }
    else if (counter[1] == 0 && iter <= adapt) {
        if ((iter >= round(adapt/5)) & (iter <= adapt)) {
            sd[1] += 1.0 * init;
            sd[2] += 1.0;
        }
    }
    if (iter == adapt) {
        sd[0] = sd[1]/sd[2];
        counter[0] = 0;
        counter[1] = 0;
    }
    return;
}

/* ************************************************************************** */

void setupproposal(long int *counter, double *sd, double init) {

    /* Set initial values to a Metropolis updater: */

    counter[0] = 0;
    counter[1] = 0;
    sd[0] = init;
    sd[1] = 0.0;
    sd[2] = 0.0;

    return;
}

/* ************************************************************************** */

void printline(double *par, int dim, FILE *file) {
    int i;
    if (dim == 1)
        fprintf(file, "%e\n", *par);
    else {
        for (i = 0; i < dim; i++)
            fprintf(file, "%e ", par[i]);
        fprintf(file, "\n");
    }
    return;
}

/* ************************************************************************** */
