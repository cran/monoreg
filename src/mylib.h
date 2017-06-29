/* ************************************************************************** */
/* File:        mylib.h                                                       */
/* Author:      Olli Saarela (olli.saarela@utoronto.ca)                       */
/* Summary:     Functions defined in file mylib.c                             */
/* ************************************************************************** */

int vidx(int o, int v, int n);

char *cvector(int length);
char **cmatrix(int rows, int cols);
double *dvector(int length);
double **dmatrix(int rows, int cols);
double **pdvector(int length);
double ***pdmatrix(int rows, int cols);
int *ivector(int length);
int **imatrix(int rows, int cols);
long int *livector(int length);
long int **limatrix(int rows, int cols);

void scrapcvector(char *p);
void scrapcmatrix(char **p, int rows);
void scrapdvector(double *p);
void scrapdmatrix(double **p, int rows);
void scrappdvector(double **p);
void scrappdmatrix(double ***p, int rows);
void scrapivector(int *p);
void scrapimatrix(int **p, int rows);
void scraplivector(long int *p);
void scraplimatrix(long int **p, int rows);

/* ************************************************************************** */
