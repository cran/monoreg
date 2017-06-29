/* ************************************************************************** */
/* File:        mylib.c                                                       */
/* Author:      Olli Saarela (olli.saarela@utoronto.ca)                       */
/* Summary:     Some functions for memory allocation etc.                     */
/* ************************************************************************** */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "mylib.h"

/* ************************************************************************** */

/* Matrix indices into vector index;
   arguments are observation, variable and number of observations: */

int vidx(int o, int v, int n) {
    return(v * n + o);
}

/* ************************************************************************** */

/* Allocate character vector: */

char *cvector(int length) {
    return(malloc(length + 1));
}

/* ************************************************************************** */

/* Allocate character matrix: */

char **cmatrix(int rows, int cols) {
    int i;
    char **p;
    p = malloc(rows * sizeof(char *));
    for(i = 0; i < rows; i++)
        p[i] = malloc(cols + 1);
    return(p);
}

/* ************************************************************************** */

/* Allocate double vector: */

double *dvector(int length) {
    return(malloc(length * sizeof(double)));
}

/* ************************************************************************** */

/* Allocate double matrix: */

double **dmatrix(int rows, int cols) {
    int i;
    double **p;
    p = malloc(rows * sizeof(double *));
    for(i = 0; i < rows; i++)
        p[i] = malloc(cols * sizeof(double));
    return(p);
}

/* ************************************************************************** */

/* Allocate pointer-to-double vector: */

double **pdvector(int length) {
    return(malloc(length * sizeof(double *)));
}

/* ************************************************************************** */

/* Allocate pointer-to-double matrix: */

double ***pdmatrix(int rows, int cols) {
    int i;
    double ***p;
    p = malloc(rows * sizeof(double *));
    for(i = 0; i < rows; i++)
        p[i] = malloc(cols * sizeof(double *));
    return(p);
}

/* ************************************************************************** */

/* Allocate integer vector: */

int *ivector(int length) {
    return(malloc(length * sizeof(int)));
}

/* ************************************************************************** */

/* Allocate integer matrix: */

int **imatrix(int rows, int cols) {
    int i;
    int **p;
    p = malloc(rows * sizeof(int *));
    for(i = 0; i < rows; i++)
        p[i] = malloc(cols * sizeof(int));
    return(p);
}

/* ************************************************************************** */

/* Allocate long integer vector: */

long int *livector(int length) {
    return(malloc(length * sizeof(long int)));
}

/* ************************************************************************** */

/* Allocate long integer matrix: */

long int **limatrix(int rows, int cols) {
    int i;
    long int **p;
    p = malloc(rows * sizeof(long int *));
    for(i = 0; i < rows; i++)
        p[i] = malloc(cols * sizeof(long int));
    return(p);
}

/* ************************************************************************** */

/* Free character vector: */

void scrapcvector(char *p) {
    free(p);
    return;
}

/* ************************************************************************** */

/* Free character matrix: */

void scrapcmatrix(char **p, int rows) {
    int i;
    for(i = 0; i < rows; i++)
        free(p[i]);
    free(p);
    return;
}

/* ************************************************************************** */

/* Free double vector: */

void scrapdvector(double *p) {
    free(p);
    return;
}

/* ************************************************************************** */

/* Free double matrix: */

void scrapdmatrix(double **p, int rows) {
    int i;
    for(i = 0; i < rows; i++)
        free(p[i]);
    free(p);
    return;
}

/* ************************************************************************** */

/* Free pointer-to-double vector: */

void scrappdvector(double **p) {
    free(p);
    return;
}

/* ************************************************************************** */

/* Free pointer-to-double matrix: */

void scrappdmatrix(double ***p, int rows) {
    int i;
    for(i = 0; i < rows; i++)
        free(p[i]);
    free(p);
    return;
}

/* ************************************************************************** */

/* Free integer vector: */

void scrapivector(int *p) {
    free(p);
    return;
}

/* ************************************************************************** */

/* Free integer matrix: */

void scrapimatrix(int **p, int rows) {
    int i;
    for(i = 0; i < rows; i++)
        free(p[i]);
    free(p);
    return;
}

/* ************************************************************************** */

/* Free long integer vector: */

void scraplivector(long int *p) {
    free(p);
    return;
}

/* ************************************************************************** */

/* Free long integer matrix: */

void scraplimatrix(long int **p, int rows) {
    int i;
    for(i = 0; i < rows; i++)
        free(p[i]);
    free(p);
    return;
}

/* ************************************************************************** */
