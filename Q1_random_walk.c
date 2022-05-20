#include<stdio.h>
#include"utility.h"
#include<math.h>
#include<stdlib.h>

#define a	572
#define m	16381



double random_walk_2D(double* xi,double* yi,int N,FILE* file){
	//file will keep record of walks
	double PI= 3.1415926;
	double x=*xi,y=*yi;
	double Xi=*xi;double Yi=*yi;
	//start walking randomly
	fprintf(file,"%lf	%lf\n",x,y);
	for(int i=1;i<=N;i++){
		x+=cos(2*PI*((1.0*mlcg_rand(a,m))/m));		
		y+=sin(2*PI*((1.0*mlcg_rand(a,m))/m));
		fprintf(file,"%lf	%lf\n",x,y);
	}

	*xi=x;*yi=y;//xi and yi will keep final x y information	
	return( sqrt(pow(x-Xi,2)+pow(y-Yi,2)) );	
}

void main(){
	//this programm is for calculations involved in q1
	//random_walk_2D() writes x,y coordinate in file
	//& returns radial displacement R 
	int seed=1234;	//randomly giving a seed for LCG random no generator
	set_seed(seed);

	FILE* file;
	FILE* file1;//for plot of Rrms vs sqrt(N)
	file1=fopen("Rrms_vs_sqrtN.txt","w");
	int K=500;//no of walks in each (steps=N)
	int N;//no of steps
	double Rrms;
	double R[K];//storing radial distance for every set of walk
	double x;
	double y;
	double del_x[K];
	double del_y[K];
	double del_X,del_Y;

	printf("\nN\tsqrt(N)\tRrms\taverage(del_X)\taverage(del_Y)\n");
N=200;
while(N<=1250){
	file=fopen("random_walk.txt","w");//this file has no use as such but  
	//function argument has file input for returning output data
	for(int i=1;i<=K;i++)
	{	x=0;y=0;
		R[i-1]=random_walk_2D(&x,&y,N,file);
		//keep x,y displacements in each set for calculation
		del_x[i-1]=x;
		del_y[i-1]=y;		
	}
	fclose(file);
	//calculate Rrms,del_X,del_Y
	del_X=0;del_Y=0;
	for(int i=1;i<=K;i++)
	{	del_X+=del_x[i-1];
		del_Y+=del_y[i-1];
		Rrms+=pow(R[i-1],2);
	}
	del_X/=K;
	del_Y/=K;
	Rrms=pow(Rrms/K,0.5);

	printf("\n%d\t%.3f\t%.3f\t%.3f\t%.3f\n",N,sqrt(1.0*N),Rrms,del_X,del_Y);
	fprintf(file1,"%lf	%lf\n",sqrt(1.0*N),Rrms);
	N=N+250;
     }

fclose(file1);	
}


//output
/*
N       sqrt(N) Rrms    average(del_X)  average(del_Y)

N       sqrt(N) Rrms    average(del_X)  average(del_Y)

200     14.142  14.419  -1.534  -0.045

450     21.213  23.050  -3.528  -0.054

700     26.458  29.841  -5.518  -0.026

950     30.822  34.727  -7.480  -0.002

1200    34.641  39.047  -9.460  0.039
*/


//comment 
/*For 200 steps 
R_rms close to sqrt(N)

but for larger values of steps there is deviation.

mlcg has a finite sequence of psudo random numbers  
for 1250 stps of each 500 walks (2*1250*500 random numbers)
probably the correlations in MLCG random numbers giving a spoil play
*/