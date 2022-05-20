#include<stdio.h>
#include<stdlib.h>

/*A4*/
double lu_dcmpsn(int n,double a[n][n],double * P,int* t);
void gauss_jordan(int n,double arr[n][n],double* b);
void partial_pivot(int n,double arr[n][n],double * P,int* t,int r);
void forward_back_sub(int n,double a[n][n],double *P,double *b);

/*A5*/
double root_newton_raphson(double(*func)(double),double x0,FILE*);
double root_false_position(double (*func)(double),double x1,double x2,FILE*);
double root_bisection(double(*func)(double),double x1,double x2,FILE *);
int brac(double (*func)(double),double* x1,double* x2);

double func(double);
double D1(double(*func)(double),double);
double D2(double(*func)(double),double);

void drive_laguerre();
double laguerre(double * arr,int n,double x0,int * iter);
double poly(double* arr,int n,double x0);
double D1_poly(double* arr,int n,double x0);
double D2_poly(double * arr,int n,double x0);

/*A6*/
double f1(double);
double f2(double);
double f3(double);
double NI_midpoint(double(*func)(double),double,double,int);
double NI_trapezoid(double(*func)(double),double,double,int);
double NI_simpson(double(*func)(double),double ,double ,int);
double monte_carlo_1D(double(*func)(double),double,double,double,double*);

/*ASSIGNMENT 7*/
void Euler_explicit(double(*func)(double,double),double x_0,double y_0,double x_n,double h,FILE * fptr);

double RK4(double(*func)(double,double,double),double x_0,double y_0,double y1_0,double x_n,double h,FILE * fptr,char ch);

double shoot(double(*func)(double,double,double),double a,double y_a,double b,double y_b,double z,double h,FILE* fptr);

/*PROJECT*/
double volume_montecarlo(int(*func)(double,double,double),double x1,double x2,double y1,double y2,double z1,double z2,double N);


/**/
double least_square_fit(FILE* file,int N,double *a0,double* a1);

/*random no genarator------*/
int mlcg_rand(int a,int m);
void set_seed(int seed);



