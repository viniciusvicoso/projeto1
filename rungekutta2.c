#include <stdio.h>
#include <math.h>

#define N 2
#define M 2
#define ti 0
#define tf 1


typedef double(*sistfunc)();


double f1(double t, double y[M])
{
	return (y[0] - y[1] + 2);
}


double f2(double t, double y[M])
{
	return (-y[0] + y[1] + 4*t); 
}


double rk(sistfunc func[], double y[M], double t, double h)	
{
	double k1[M], k2[M], k3[M], k4[M], yp[M];
	int i;
	
	
	for(i = 0; i < M; i++)
		k1[i] = h*func[i](t, y);
	
	for( i = 0; i < M; i++)
		yp[i] = y[i] + (1/2.0)*k1[i];
	
	for(i = 0; i < M;i++)
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
	
	
	return y[1];
}	
	
void main(int arg, char **arqv)
{
	FILE *out;
	double u0[M] = {-1,0}, u[M], h, t=ti, NP = 10.0;
	int i;
	sistfunc equacoes[M] = {f1,f2};
	
	
	out = fopen(arqv[1], "w");
	
	h = (tf-ti) / NP; 
	
	
	for(i = 0; i < M; i++)
		u[i] = u0[i];
			
			
	do
	{
		fprintf(out, "%lf\t%lf\n", t, u[1]);
		
		u[1] = rk(equacoes,u,t,h);
		
		t += h;
	
	}while(t <= tf);
	
	fprintf(out, "%lf\t%lf\n", t, u[1]);
}
