#include <stdio.h>
#include <math.h>

#define M 3
#define ti 0
#define tf 3


typedef double(*sistfunc)();


double f1(double t, double y[M])
{
	return(y[2]);
}


double f2(double t, double y[M])
{
	return(y[1]);
}


double f3(double t, double y[M]) 
{
	return(exp(t) - 2*y[2] + y[1] + 2*y[0]);
}


double rk(sistfunc func[], double y[M], double t, double h)	
{
	double k1[M], k2[M], k3[M], k4[M], yp[M];
	int i;
	
	
	for(i = 0; i < M; i++)	
		k1[i] = h*func[i](t, y);
		
	for(i = 0; i < M; i++)
		yp[i] = y[i] + (1/2.0)*k1[i];
		
	for(i = 0; i < M; i++)
		k2[i] = h*func[i](t + 0.5*h, yp);
		
	for(i = 0; i < M; i++)
		yp[i] = y[i] + (1/2.0)*k2[i];
		
	for(i = 0; i < M; i++)
		k3[i] = h*func[i](t + 0.5*h, yp);
		
	for(i = 0; i < M; i++)
		yp[i] = y[i] + k3[i];
		
	for(i = 0; i < M; i++)
		k4[i] = h*func[i](t + h, yp);
	
	for(i = 0; i < M; i++)
		y[i] += (1/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
	
	
	return (y[0]);
}	

	
void main(int arg, char **arqv)
{
	FILE *in;
	double y0[M] = {1,2,0}, y[M], h, t = ti, N = 30.0;
	int i;
	sistfunc equacoes[M] = {f1,f2,f3};
	
	
	in = fopen(arqv[1], "w");
	
	h = (tf-ti)/N;
	
	
	for(i = 0; i < M; i++)
		y[i] = y0[i];
	
	do
	{
		fprintf(in, "%lf\t%lf\n", t, y[0]);
		
		y[0] = rk(equacoes,y,t,h);
		
		t += h;
	
	}while(t <= tf);
	
	fprintf(in, "%lf\t%lf\n", t, y[0]);
}
