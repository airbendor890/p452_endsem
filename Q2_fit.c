#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"utility.h"

//Chi	SQUARE FITTING OF DATA
//for least squares take sigma =1


/*fit invloving these choice of basis functions(legendre polynomial)*/
//y= [a_0]*P_0(x)+[a_1]*P_2(x)+[a_2]*P_4(x)+[a_3]*P_6(x)+
//double P_6(double z){return  (1.0/16) * ( 231*pow(z,6)-315*pow(z,4)+105*pow(z,2)-5) ; }
double P_4(double z){return (1.0/8)*(35*z*z*z*z-30*z*z+3);}
double P_3(double z){return 0.5*(5*z*z*z-3*z);}
double P_2(double z){return 0.5*(3*z*z-1);}
double P_1(double z){return z;};
double P_0(double z){return 1;}
int main(){
	
/*number of basis functions M,number of data points N,Matrix A_N*M ,
x[N]:x points,y[N]:y points
*/

	
	int N=26;					//Number of Data points
	int M=5;					//Number of Parameters
	double a[M];				//parameters
	double A[N][M];				//Design Matrix
	double x[N];
	double y[N];
	double sig[N];				//normal variance at each y_i
	double chisq;
	double (*fn_ptr_arr[])(double)={P_0,P_1,P_2,P_3,P_4};//array of function pointers

	/*read data points from input file*/
	
	FILE* file1;
	file1=fopen("esem4fit.txt","r");

	for(int i=0;i<N;i++){
		fscanf(file1,"%lf",&x[i]);
		fscanf(file1,"%lf",&y[i]);
		sig[i]=1;
		

	}
	
	fclose(file1);


	//Make thr design matrix of fitting problem
	for(int i=0;i<N;i++){
	for(int j=0;j<M;j++)
		A[i][j]=(*fn_ptr_arr[j])(x[i])/sig[i];
	}




	/*matrix [alpha]=A^t*A,M*M;matrix C=inerse(alpha);C is the covariance matrix*/
	double alpha[M][M];
	//[alpha]=A^t*A
	for(int i=0;i<M;i++){
		for(int j=0;j<M;j++){
			double sum=0;
			for(int k=0;k<N;k++)
				sum+=A[k][i]*A[k][j];
			alpha[i][j]=sum;
		}

	}

	//covariance matrix C=inverse([alpha])
	//Using Gauss jordan method
	long double C[M][M];

	for(int i=0;i<M;i++){
		double e[M];
		for(int j=0;j<M;j++){
			if(j==i)	e[j]=1;
			else		e[j]=0;
		}
		//Work with copy of alpha matrix since gauss jordan elemination
		//will change the matrix
		double copy_alpha[M][M];
		
		for(int j=0;j<M;j++)
		for(int k=0;k<M;k++)
			copy_alpha[j][k]=alpha[j][k];
	
		//call gauss jordan
		gauss_jordan(M,copy_alpha,e);

		for(int k=0;k<M;k++)
		C[k][i]=e[k];	
	}


	//check
	/*
	printf("\n\n");
	for(int i=0;i<M;i++){
		for(int j=0;j<M;j++){
			printf("%lf\t",C[i][j]);
		}
		printf("\n");	
	}
	*/
	
	
	double b[M];
	for(int i=0;i<M;i++){
		double sum=0;
		for(int k=0;k<N;k++)
			sum+=A[k][i]*(y[k]/sig[k]);
		b[i]=sum;
	}


	//evaluate parameter vector
	for(int i=0;i<M;i++){
		double sum=0;
		for(int k=0;k<M;k++)
			sum+=C[i][k]*b[k];
		a[i]=sum;
	}

	printf("\nparameter\tvalue\terror\n");
	for(int i=0;i<M;i++){
		printf("%d\t%le\t%le\n",i,a[i],sqrt(C[i][i]));	
	
	}
	
	printf("\nCovariance Matrix::\n");
	for(int i=0;i<M;i++){
		printf("\n");
		for(int j=0;j<M;j++){
			printf("%Le\t",C[i][j]);
		}
		printf("\n");
	}	
	
	//Evaluate Chi Square
	
	chisq=0.0;
	for(int i=0;i<N;i++){
		double temp=0.0;
		for(int j=0;j<M;j++){
			temp+=a[j]*(*fn_ptr_arr[j])(x[i]);
		}
		chisq+=pow((y[i]-temp)/sig[i],2);
	}
	
	printf("\nchisq/dof = %lf\n",chisq/(N-M));
	
	return 0;
}


//output

/*

parameter       value   error
0       6.965780e-02    1.976653e-01
1       3.624020e-03    3.313630e-01
2       -1.208258e-02   4.181216e-01
3       1.142622e-02    4.661041e-01
4       1.104924e-01    5.143079e-01

Covariance Matrix::

3.907156e-02    -4.482551e-18   -4.841037e-03   7.011137e-19    -9.553491e-03

-4.482551e-18   1.098014e-01    9.826637e-18    -2.537516e-02   -2.988966e-18

-4.841037e-03   9.826637e-18    1.748257e-01    -8.691486e-18   -4.937569e-02

7.011137e-19    -2.537516e-02   -8.691486e-18   2.172530e-01    1.426517e-17

-9.553491e-03   -2.988966e-18   -4.937569e-02   1.426517e-17    2.645126e-01

chisq/dof = 0.000134
*/