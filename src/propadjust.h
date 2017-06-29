/* ************************************************************************** */
/* File:        propadjust.h                                                  */
/* Author:      Olli Saarela (olli.saarela@utoronto.ca)                       */
/* Summary:     Functions defined in the file propadjust.c                    */
/* ************************************************************************** */

double sdadj(double sd, double rate);
void adjustproposal(long int *counter, double *sd, int iter, int adapt, double init);
void setupproposal(long int *counter, double *sd, double init);
void printline(double *par, int dim, FILE *file);
