#include<stdio.h>
#include<math.h>

double func(double x){
	double PI= 3.1415926;
	return 20*fabs(sin(PI*x));
}

//1-D diffusion PDE :Explicit
void diffusion_1D_explicit(double(*func)(double),double dx,double dt,double x1,double x2,int N_T,double * V){
//space domain:[a,b]
//discrete element:space=>dx,time=>dt
int N=(x2-x1)/dx;
double alpha=dt/(dx*dx);
double V_0[N+1];

for(int i=0;i<=N;i++){
	V_0[i]=(*func)(x1+dx*i);
	*(V+i)=V_0[i];
}
for(int i=0;i<N_T;i++){
	for(int j=1;j<N;j++){
		*(V+j)=(1-2*alpha)*V_0[j]+alpha*(V_0[j+1]+V_0[j-1]);
		
	}
	
	for(int j=1;j<N;j++){
		V_0[j]=*(V+j);
	}	
}
//for(int j=0;j<=N;j++)printf("%lf\n",*(V+j));
}

void main(){
	
	double dx=0.1;
	double dt=0.0008;
	double x1=0;
	double x2=2;
	
	//time t=0
	int  n_x=((x2-x1)/dx);
	printf("n_x%d",n_x);
	double V[n_x+1];
	diffusion_1D_explicit(&func,dx,dt,x1,x2,0,V);
	
	FILE* file=fopen("metal_rod_00.txt","w");
	for(int i=0;i<=n_x;i++){
		fprintf(file,"%lf\t%lf\n",i*dx,V[i]);
	}
	fclose(file);
	
	//time t=10*dt
	diffusion_1D_explicit(&func,dx,dt,x1,x2,10,V);
	
	file=fopen("metal_rod_10.txt","w");
	for(int i=0;i<=n_x;i++){
		fprintf(file,"%lf\t%lf\n",i*dx,V[i]);
	}
	fclose(file);

	//time t=20*dt
	diffusion_1D_explicit(&func,dx,dt,x1,x2,20,V);
	
	file=fopen("metal_rod_20.txt","w");
	for(int i=0;i<=n_x;i++){
		fprintf(file,"%lf\t%lf\n",i*dx,V[i]);
	}
	fclose(file);	

	//time t=50*dt
	diffusion_1D_explicit(&func,dx,dt,x1,x2,50,V);
	
	file=fopen("metal_rod_50.txt","w");
	for(int i=0;i<=n_x;i++){
		fprintf(file,"%lf\t%lf\n",i*dx,V[i]);
	}
	fclose(file);
	
	//time t=100*dt
	diffusion_1D_explicit(&func,dx,dt,x1,x2,100,V);
	
	file=fopen("metal_rod_100.txt","w");
	for(int i=0;i<=n_x;i++){
		fprintf(file,"%lf\t%lf\n",i*dx,V[i]);
	}
	fclose(file);
	
	
	//time t=200*dt
	diffusion_1D_explicit(&func,dx,dt,x1,x2,200,V);
	
	file=fopen("metal_rod_200.txt","w");
	for(int i=0;i<=n_x;i++){
		fprintf(file,"%lf\t%lf\n",i*dx,V[i]);
	}
	fclose(file);
	
	//time t=500*dt
	diffusion_1D_explicit(&func,dx,dt,x1,x2,500,V);
	
	file=fopen("metal_rod_500.txt","w");
	for(int i=0;i<=n_x;i++){
		fprintf(file,"%lf\t%lf\n",i*dx,V[i]);
	}
	fclose(file);	
}