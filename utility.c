#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"utility.h"



/*--------------------------------A4-----------------------------------------*/

//-----------------------------------------------------------------------
// function for partial pivoting
void partial_pivot(int n,double arr[n][n],double* P,int* t,int r){
/*when using for gauss jordan P is vector b*/
/*when using for LU dcmpsn ,P is for storing permutaion info(vector:1,2...,n)*/
/* t will hold no of permutation that wil occur .useful in det calculation*/

		int r1=r+1;
		while(arr[r][r]==0 && r1<n){
			if(abs(arr[r1][r])>abs(arr[r][r])){
				double temp;
				for(int c=0;c<n;c++){
					temp=arr[r1][c];
					arr[r1][c]=arr[r][c];
					arr[r][c]=temp;
				}

				int dum;
				dum=P[r1];
				P[r1]=P[r];
				P[r]=dum;

				*t=*t+1;	//row exchange occuured
			}
			r1++;
		}



}
/*---------------------------------------------------------------------------*/
/*gauss jordan elimination function*/

void gauss_jordan(int n,double arr[n][n],double* b){
	int t=0;//for storing no of row exchange in partial pivot
	for(int r=0;r<n;r++){
		if(abs(arr[r][r])<1.0e-12)
		partial_pivot(n,arr,b,&t,r);
		double pivot=arr[r][r];
		for(int c=r;c<n;c++)
		arr[r][c]=arr[r][c]/pivot;
		b[r]=b[r]/pivot;

		for(int r1=0;r1<n;r1++){
			if(r1==r || arr[r1][r]==0)continue;
			else{
				double factor=arr[r1][r];
				for(int c=r;c<n;c++)
				arr[r1][c]=arr[r1][c]-factor*arr[r][c];
				b[r1]=b[r1]-factor*b[r];
													
			}
		}

	}
}

//-----------------------------------------------------------------------
//lu dcmpsn function
double lu_dcmpsn(int n,double a[n][n],double *P,int* t){
	/* P will hold the permutation infromation that will occur during 
	partial pivoting*/
	/*t will hold no of permutation that wil occur .useful in det calculation*/

	double sum=0.0;
	double det=1.0;
	if(a[0][0]==0) partial_pivot(n,a,P,t,0);

	for(int j=0;j<n;j++){	/*loop over columns of crouts method*/
		for(int i=0;i<=j;i++){
			sum=a[i][j];
			for(int k=0;k<i;k++) sum=sum-a[i][k]*a[k][j];
			a[i][j]=sum;
			
		}
		

		for(int i=j+1;i<n;i++){
			sum=a[i][j];
			for(int k=0;k<j;k++) sum=sum-a[i][k]*a[k][j];
			//still if diag element is 0 its a singular matrix
			//return 0  
			if(a[j][j]==0.0){	
				printf("cant divide by 0! singular matrix !!\n LU decomposition failur!!\n");
				return 0;
			}
			else
			a[i][j]=sum/a[j][j];
		}

	}

	//its not singular ..return its determinant
	for(int i=0;i<n;i++)det=det*a[i][i];
	det=det*pow(-1,*t);

	return det;

}	

//-------------------------------------------------------------------

void forward_back_sub(int n,double a[n][n],double *P,double * b){
	/* first make the permutation matrix from P*/
	int perm[n][n];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(j==P[i])
				perm[i][j]=1;
			else
				perm[i][j]=0;

		}
	}

	/* now permute the vector b*/
	double c[n];
	for(int i=0;i<n;i++){
		double sum=0.0;
		for(int j=0;j<n;j++)
			sum+=perm[i][j]*b[j];
		c[i]=sum;
	}
	for(int i=0;i<n;i++)b[i]=c[i];

	
	double sum=0.0;
	/*forward substitution Ly=b */
	for(int i=0;i<n;i++){
		sum=b[i];
		for(int j=0;j<i;j++){
			sum=sum-a[i][j]*b[j];
		}
		b[i]=sum;
	}

	/* back subtitution Ux=y */

	for(int i=n-1;i>=0;i--){
		sum=b[i];
		for(int j=i+1;j<n;j++){
			sum=sum-a[i][j]*b[j];
		}

		b[i]=sum/a[i][i];
	}

}

//------------------------------------------------------------


/*----------------------A5-------------------------------------------------------*/


//--------------------------------------------------------------------------------
/*root bracketing function*/
int brac(double (*func)(double),double* x1,double* x2){
	double f1,f2;
	double FACTOR=1.5;
	
	//if(*x1=*x2)break;/*bad range*/
	f1=(*func)(*x1);
	f2=(*func)(*x2);
	for(int j=1;j<=50;j++){
		if(f1*f2<0.0) return 1;	/*change of sign...bracketed*/		
		if(fabs(f1)<fabs(f2))
			f1=(*func)(*x1+=FACTOR*(*x1-*x2));
		else
			f2=(*func)(*x2+=FACTOR*(*x2-*x1));
	}

	return 0;	/*not bracketed even after 50 trys*/
}


//---------------------------------------------------------------------------------
/*root bisection method function*/
double root_bisection(double (*func)(double),double x1,double x2,FILE* fptr)
{	/*assuming input is such that x2>x1*/
	/* fptr will write the iterations to file*/
	/*it must be passed with file opened*/
	double xmid;	
	int JMAX=200;
	double epsilon=0.000001;

	for(int j=1;j<=JMAX;j++){
		xmid=0.5*(x1+x2);

		if((*func)(x1)*(*func)(xmid)<0)
			x2=xmid;
		else			
			x1=xmid;
		
		printf("\n%d	%lf	%lf	%lf\n",j,x1,x2,fabs(x1-x2));
		fprintf(fptr,"%d	%lf\n",j,fabs(x1-x2));
		/*if converged to desirable accuracy*/
		if(fabs(x1-x2)<epsilon){fclose(fptr);return xmid;}

	}

}

//------------------------------------------------------------------------------------
/*false position root finding method function*/
double root_false_position(double (*func)(double),double x1,double x2,FILE * fptr)
{
	/*assumption is that x1 and x2 are passsed after bracketing*/
	/*fptr will write the data to opened file*/
	/*assumption is that file is already open before passing*/
	double epsilon=0.0001;
	int JMAX=200;
	double c_k,c_k1;

	c_k1=x2-((x2-x1)*(*func)(x2))/((*func)(x2)-(*func)(x1));
	if((*func)(x1)*(*func)(c_k1)<0)
		x2=c_k1;
	else
		x1=c_k1;


	for(int j=1;j<=JMAX;j++){
		c_k=c_k1;
		c_k1=x2-((x2-x1)*(*func)(x2))/((*func)(x2)-(*func)(x1));
		
		if((*func)(x1)*(*func)(c_k1	)<0)
			x2=c_k1;
		else
			x1=c_k1;

		printf("\n%d	%lf	%lf	%lf\n",j,c_k,c_k1,fabs(c_k1-c_k));
		fprintf(fptr,"%d	%lf\n",j,fabs(c_k1-c_k));

		if(fabs(c_k1-c_k)<epsilon) return c_k1;
	}

}

//--------------------------------------------------------------------------------------
/*newton raphson root finding method function*/
double root_newton_raphson(double(*func)(double),double x0,FILE * fptr){
	/*x0 is the gussed root*/
	/*fptr will write convergence data to file*/
	/*file pointer must be passed already linked to a file*/
	double epsilon=0.000001;
	int JMAX=40;
	double x_n,x_n1;

	x_n1=x0;
	for(int j=1;j<=JMAX;j++){
		x_n=x_n1;
		x_n1=x_n-((*func)(x_n)/D1(func,x_n));	/*first derivative*/
		fprintf(fptr,"%d	%lf\n",j,fabs(x_n1-x_n));
		printf("\n%d\t%lf\t%lf\t%lf\n",j,x_n,x_n1,fabs(x_n1-x_n));
		if(fabs(x_n1-x_n)<epsilon) return x_n;
	}

}

//-------------------------------------------------------------------------------------
/*1st derivative for one variable function*/
double D1(double(*func)(double),double x0){
	double h=0.001;
	
	double d=((*func)(x0+h)-(*func)(x0-h))/(2*h);
	return d;
}

//------------------------------------------------------------------------------------

/*2nd derivative for  one variavble function*/

double D2(double(*func)(double),double x0){
	double h=0.0001;
	double d=((*func)(x0+h)+(*func)(x0-h)-2*(*func)(x0))/(h*h);

	return d;

}

//------------------------------------------------------------------------------------
void drive_laguerre(){
	int deg;
	FILE* file;
	printf("\nLaguerre and synthetic division method for real roots of polynomial of degree n\n");
	printf("\nEnter degree here and store coefficients of polynomial in a text file(poly.txt) in descending power:\n");
	scanf("%d",&deg);
	double arr[deg+1];
	double roots[deg];/*store roots*/
	int iter;
	file=fopen("poly.txt","r");
	for(int i=0;i<=deg;i++)
		fscanf(file,"%lf",&(arr[i]));
	fclose(file);	

	printf("\nroot		iterations taken\n");
	for(int j=0;j<deg;j++){
		double x0=5.00;
		roots[j]=laguerre(arr,deg-j,x0,&iter);
		printf("\n%lf		%d\n",roots[j],iter);
	
		/*perform deflation*/
		double r=0.0;
		for(int jj=0;jj<=deg-j;jj++){
			arr[jj]+=r*(roots[j]);
			r=arr[jj];

		}
		/*end deflation*/		
	}


}

//------------------------------------------------------------------------------------
/*function: laguerre method to find a real root of a polynomial*/

double laguerre(double * arr,int n,double x0,int * iter){
	/*array input n degree polynomial,guess root x0*/
	/*no of iteration taken in *iter */
	int MAXIT=50;
	double epsilon=1.0e-6;
	if(fabs(poly(arr,n,x0))<epsilon){*iter=0;return x0;}
	double G,H,a_k,a_k1;
	
	a_k1=x0;
	for(int j=1;j<MAXIT;j++){
		double a,d1,d2;
		a_k=a_k1;

		G=D1_poly(arr,n,a_k)/poly(arr,n,a_k);
		H=pow(G,2)-(D2_poly(arr,n,a_k)/poly(arr,n,a_k));
		d1=G+sqrt((n-1)*(n*H-pow(G,2)));
		d2=G-sqrt((n-1)*(n*H-pow(G,2)));
	
		if(fabs(d1)>fabs(d2))
			a=n/d1;
		else
			a=n/d2;
				
		a_k1=a_k-a;
		//printf("\n%lf	%lf	%lf\n",a,d1,d2);
		if(fabs(a_k1-a_k)<epsilon){*iter=j; return a_k;}
		if(fabs(poly(arr,n,a_k1))<epsilon){*iter=j;return a_k1;}

	}
	
	*iter=MAXIT;/*failed to find root within MAXIT*/

}

//------------------------------------------------------------------------------------
/*function:polynomial value at x0(array input)*/
double poly(double* arr,int n,double x0){
	double value=0.0;
	for(int i=0;i<=n;i++)
		value+=arr[i]*pow(x0,n-i);
	return value;
}

//------------------------------------------------------------------------------------
/*function: first derivarive of a polynomial function At x0(array input)*/
double D1_poly(double* arr,int n,double x0){
	double value=0.0;
	for(int i=0;i<n;i++)
		value+=arr[i]*(n-i)*pow(x0,n-i-1);
//	printf("\n%lf\n",value);
	return value;
}

//------------------------------------------------------------------------------------
/*function: second derivative of a polynomial function at x0(array input)*/
double D2_poly(double * arr,int n,double x0){
	double value=0.0;
	for(int i=0;i<n-1;i++)
		value+=arr[i]*(n-i)*(n-i-1)*pow(x0,n-i-2);

	return value;
}

//-----------------------------------------------------------------------------------
//-------------------------A6--------------------------------------------------------
/*NUmerical integration midpoint*/

double NI_midpoint(double(*func)(double),double a,double b,int N){
	double sum,h;
	h=(b-a)/N;
	sum=0.0;
	double x=a;
	for(int j=1;j<=N;j++){
		sum+=((*func)(x+(j-1)*h)+(*func)(x+j*h))/2;
	}
	return h*sum;
}
/*---------------------------------------------------------------------------------*/
/*numerical integration	trapezoid*/
double NI_trapezoid(double(*func)(double),double a,double b,int N){
	double sum,h;
	h=(b-a)/N;
	sum=0.0;
	for(int j=0;j<=N;j++){
		if(j==0 || j==N)
			sum+=(*func)(a+j*h);
		else
			sum+=2*(*func)(a+j*h);
	}
	return (h*sum)/2.0;
}
/*-------------------------------------------------------------------------------*/
/*numerical integration simpson*/
double NI_simpson(double(*func)(double),double a,double b,int N){
	if(N%2!=0)N=N+1;
	double sum,h;
	h=(b-a)/(N);/*dividing into even no of interval*/
	sum=0.0;
	for(int j=0;j<=N;j++){
		if(j==0 || j==N)
			sum+=(*func)(a+j*(h));
		else if(j%2==0)
			sum+=2.0*(*func)(a+j*(h));
		else
			sum+=4.0*(*func)(a+j*(h));
	}
	return (h*sum)/3.0;
}
/*------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------*/
/*Numerical integration monte-carlo method */
double monte_carlo_1D(double(*func)(double),double a,double b,double N,double* s){
						/*return sigma in s*/
	double F_N=0;
	double sum1=0;
	double sum2=0;
	double x;
	
	time_t t;

	srand((unsigned) time(&t));

	for(int i=1;i<=N;i++){
		x=a+(b-a)*((1.0*rand())/RAND_MAX);
		sum1+=(*func)(x);
		sum2+=pow((*func)(x),2);
	}
	F_N=((b-a)/N)*sum1;
	*s=sqrt((sum2/N)-pow(sum1/N,2));
	return F_N;
}
/*----------------------------------------------------------------------------*/
/*Assignment 7:numerical methods solving of differentila eqn*/
/*Eulers explicit method:*/

void Euler_explicit(double(*func)(double,double),double x_0,double y_0,double x_n,double h,FILE * fptr)
/*x_0,y_0 is initial value,x_n-x_0 is the interval of plot, h is step size,fptr is pointer to output file */
{	
	int k=(x_n-x_0)/h;
	double x=x_0;
	double y=y_0;
	fprintf(fptr,"%lf	%lf\n",x,y);

	for(int i=1;i<=k;i++){
		y=y+h*(*func)(y,x);
		x=x+h;
		fprintf(fptr,"%lf	%lf\n",x,y);
	}

}
/*--------------------------------------------------------------------------*/
/*Runge kutta method RK4*/
double RK4(double(*func)(double,double,double),double x_0,double y_0,double y1_0,double x_n,double h,FILE * fptr,char ch)
/*this RK4 routine is for solving 2nd order ode with initial values (or twocoupled 1st order ode)*/
/*x_0,y_0,y1_0 are initial value,x_n-x_0 is the interval of plot, h is step size,fptr is pointer to output file */
/* func:=f(dy/dx,y,x)*/
/*pass ch='r' if output is to be written in file*/
{	printf("\nrk4");
	int k=(x_n-x_0)/h;
	double x=x_0;
	double y=y_0;/*y(x_0)=y_0*/
	double y1=y1_0;/*y'(x_0)=y1_0*/
	fprintf(fptr,"%lf\t%lf\t%lf\n",x,y,y1);
	double k1_y,k2_y,k3_y,k4_y;
	double k1_y1,k2_y1,k3_y1,k4_y1;
	for(int i=1;i<=k;i++)
	{
		k1_y=h*y1;
		k1_y1=h*(*func)(y1,y,x);

		k2_y=h*(y1+k1_y1/2);
		k2_y1=h*(*func)(y1+k1_y1/2,y+k1_y/2,x+h/2);

		k3_y=h*(y1+k2_y1/2);
		k3_y1=h*(*func)(y1+k2_y1/2,y+k2_y/2,x+h/2);

		k4_y=h*(y1+k3_y1);
		k4_y1=h*(*func)(y1+k3_y1,y+k3_y,x);

		y=y+(k1_y+2*k2_y+2*k3_y+k4_y)/6;
		y1=y1+(k1_y1+2*k2_y1+2*k3_y1+k4_y1)/6;
		x=x+h;
	if(ch=='r')/*write if char ch is r*/
	fprintf(fptr,"%lf\t%lf\t%lf\n",x,y,y1);

	}
	return y;/*return y at x_n(used in solving boundary value 		 problem)*/		
}

/*Shooting method for boundary value problem(2nd order)*/
double shoot(double(*func)(double,double,double),double a,double y_a,double b,double y_b,double z,double h,FILE* fptr)
/*this shoot routine is for solving 2nd roder ode with boundary values*/
/*uses RK4 routine*/
/*y_a,y_b are boundary values at x=a and x=b respectively,h is step size,fptr is pointer to output file*/
/* z :guessed slope at x=a */
{
	printf("\nshoot");
	double z_h,z_l,beta,beta_z_h,beta_z_l;
	beta=RK4(func,a,y_a,z,b,h,fptr,'n');/*n :dont write in file*/
	beta_z_l=beta;
	beta_z_h=beta;
	if(fabs(beta-y_b)<0.001){

		return RK4(func,a,y_a,z,b,h,fptr,'r');/*r:write in file*/
		
	}
	else {
		printf("\nstill shooting");
		if(beta>y_b){
			z_h=z;
			while(beta>y_b){
				z=z-1.5;
				beta=RK4(func,a,y_a,z,b,h,fptr,'n');
			}
			z_l=z;
			beta_z_l=beta;
			/*lagrange interpolation*/
			z=z_l+(z_h-z_l)*(y_b-beta_z_l)/(beta_z_h-beta_z_l);
			z=(z_l*(beta_z_h-y_b)+z_h*(y_b-beta_z_l))/(beta_z_h-beta_z_l);
			shoot(func,a,y_a,b,y_b,z,h,fptr);
		}
		else{
			z_l=z;
			while(beta<y_b){
				z=z+1.5;
				beta=RK4(func,a,y_a,z,b,h,fptr,'n');
			}
			z_h=z;
			beta_z_h=beta;
			/*lagrange interpolation*/
			z=z_l+(z_h-z_l)*(y_b-beta_z_l)/(beta_z_h-beta_z_l);

			shoot(func,a,y_a,b,y_b,z,h,fptr);

		}
	}
		
}

/*-----------------------------------------------------------------------*/
//Project`````````````````````````````````````````````//
double volume_montecarlo(int(*func)(double,double,double),double x1,double x2,double y1,double y2,double z1,double z2,double N)
/*this function will evaluate volume of intrest using monte carlo method*/
/*func will return 0:(x,y,z) lies outside volume of intrest
		   1:(x,y,z) lies inside volume of intrest
*/
{
	double F_N=0;
	double sum=0;
	double x,y,z;
	FILE* file;
	file=fopen("volume_montecarlo_data.txt","w");	
	/*points that occur inside the volume of intrest will be collected to file for plotting*/	

	for(int i=1;i<=N;i++){

		/*random point in 3-D*/
		x=x1+(x2-x1)*((1.0*rand())/RAND_MAX);
		y=y1+(y2-y1)*((1.0*rand())/RAND_MAX);
		z=z1+(z2-z1)*((1.0*rand())/RAND_MAX);
		/*if the point is inside desired volume it will be counted*/
		if((*func)(x,y,z)==1)
		{	sum+=1;
			fprintf(file,"%lf\t%lf\t%lf\n",x,y,z);
		}
	}
	F_N=(((x2-x1)*(y2-y1)*(z2-z1))/N)*sum;
	return F_N;
}

/*double random_walk_2D(double* xi,double* yi,int N,FILE* file){
	//file will keep record of walks
	double PI= 3.1415926;
	double x=*xi,y=*yi;
	double Xi=*xi;double Yi=*yi;
	//start walking randomly
	fprintf(file,"%lf	%lf\n",x,y);
	for(int i=1;i<=N;i++){
		x+=cos(2*PI*((1.0*rand())/RAND_MAX));		
		y+=sin(2*PI*((1.0*rand())/RAND_MAX));
		fprintf(file,"%lf	%lf\n",x,y);
	}

	*xi=x;*yi=y;//xi and yi will keep final x y information	
	return( sqrt(pow(x-Xi,2)+pow(y-Yi,2)) );	
}


*/
/*---------------------End sem---------------*/
double least_square_fit(FILE* file,int N,double *a0,double* a1)
//file:data file on which fitting will be done
//N:no of data points
//parameters:a0,a1

{	double x[N];
	double y[N];
	for(int i=0;i<N;i++){
		fscanf(file,"%lf\t%lf",&(x[i]),&(y[i]));
	}	
	double S=N;
	double Sx,Sy,Sxx,Syy,Sxy;/*all have usual meaning*/
	double r;	/*pearsons r*/
	double delta;	
	Sx=0;Sy=0;Sxx=0;Syy=0;Sxy=0;
	
	for(int i=0;i<N;i++){
		Sx=Sx+x[i];
	}

	for(int i=0;i<N;i++){
		Sy+=y[i];
	}
	
	for(int i=0;i<N;i++){
		Sxx+=pow(x[i],2);
	}
	
	for(int i=0;i<N;i++){
		Syy+=pow(y[i],2);
	}
	

	for(int i=0;i<N;i++){
		Sxy+=x[i]*y[i];
	}
	delta=S*Sxx-pow(Sx,2);
	*a0=(Sxx*Sy-Sx*Sxy)/delta;
	*a1=(S*Sxy-Sx*Sy)/delta;			
	
	r=sqrt(pow(Sxy,2)/(Sxx*Syy));
	/*return pearsons r*/	
	return r;
}


//Random number generator

//multiplicative linear congruential generator
int prev_rand;
void set_seed(int seed){
	prev_rand=seed;
}
int mlcg_rand(int a,int m){
	int next_rand=(a*prev_rand)%m;
	prev_rand=next_rand;
	return next_rand;
}

