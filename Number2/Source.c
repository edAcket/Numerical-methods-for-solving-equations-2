#include <stdio.h>     
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <process.h>
#include <time.h>

double hx,hy,tau,*A,*B,*C,*F;
int i, j, NX = 50, NY = 50, NT = 50, NN = 1000;
const double PI = 3.141592653589793;

double** progonX(int k,double **U0)
 {
 double g1 = tau / (hx*hx);
 double g2 = tau / (hy*hy);
 double *d, *sigma,*U,**U1;
 U1 = (double**)malloc(NX * sizeof(double*));
 for (i = 0; i < NX; i++)  // цикл по строкам
	 U1[i] = (double*)malloc(NY * sizeof(double));
 A = (double*)malloc(NN * sizeof(double));
 B = (double*)malloc(NN * sizeof(double));
 C = (double*)malloc(NN * sizeof(double));
 U = (double*)malloc(NN * sizeof(double));
 F = (double*)malloc(NN * sizeof(double));
 d = (double*)malloc(NN * sizeof(double));
 sigma = (double*)malloc(NN * sizeof(double));
 for (int i = 0; i < NX; i++)
 {
	 U1[i][0] = U0[i][0];
	 U1[i][NY-1] = U0[i][NY-1];
 }
 	 for (int j = 1; j < NY - 1; j++)
		{
		 for (int i = 1; i < NX - 1; i++)
			 {
				 A[i] = 0.5 * g1;
				 C[i] = 1 + g1;
				 B[i] = 0.5 * g1;
				 F[i] = 0.5*g2*U0[i][j - 1] + (1 - g2) * U0[i][j] + 0.5*g2*U0[i][j + 1];
			 }
			 sigma[0] = 0;
			 d[0] = 0;
			  for (int i = 1; i < NX - 1; i++)
					 {
						 d[i] = B[i] / (C[i] - d[i-1] * A[i]);
						 sigma[i] = (A[i] * sigma[i-1] + F[i]) / (C[i] - d[i-1]* A[i]);
					 }
			  U[NX-1] = sigma[NX-2] / (1 - d[NX-2]);
			  for (int i = NX - 2; i >= 0; i--)
				  {
				  U[i] = d[i]* U [i + 1] + sigma[i];
				  }
		 for (int i = 0; i < NX; i++)
			 U1[i][j] = U[i];
		 }
	 return (U1);
	 free(U1), free(U), free(A); free(B); free(C); free(F); free(d); free(sigma);
 }
double** progonY(int k, double **U0)
{
	double g1 = tau / (hx*hx);
	double g2 = tau / (hy*hy);
	double *d, *sigma, *U, **U1;
	U1 = (double**)malloc(NX * sizeof(double*));
	for (i = 0; i < NX; i++)  // цикл по строкам
		U1[i] = (double*)malloc(NY * sizeof(double));
	A = (double*)malloc(NN * sizeof(double));
	B = (double*)malloc(NN * sizeof(double));
	C = (double*)malloc(NN * sizeof(double));
	U = (double*)malloc(NN * sizeof(double));
	F = (double*)malloc(NN * sizeof(double));
	d = (double*)malloc(NN * sizeof(double));
	sigma = (double*)malloc(NN * sizeof(double));
	for (int j = 0; j < NY; j++)
	{
		U1[0][j] = U0[0][j];
		U1[NX-1][j] = U0[NX-1][j];
	}
	for (int i = 1; i < NX - 1; i++)
	{
		for (int j = 1; j < NY - 1; j++)
		{
			A[j] = 0.5 * g2;
			C[j] = 1 + g2;
			B[j] = 0.5 * g2;
			F[j] = 0.5*g1*U0[i-1][j] + (1 - g1) * U0[i][j] + 0.5*g1*U0[i+1][j];
		}
		sigma[0] = 0;
		d[0] = 0;
		for (int j = 1 ; j < NY - 1; j++)
		{
			d[j] = B[j] / (C[j] - d[j-1] * A[j]);
			sigma[j] = (A[j] * sigma[j-1] + F[j]) / (C[j] - d[j-1] * A[j]);
		}
		U[NY-1] = 0;
		for (int j = NY - 2; j >= 0; j--)
		{
			U[j] = d[j] * U[j + 1] + sigma[j];
		}
		for (int j = 0; j < NY; j++)
			U1[i][j] = U[j];
	}
	return (U1);
	free(U1), free(U), free(A); free(B); free(C); free(F); free(d); free(sigma);
}
int main(void)
{

	FILE *fp1;
	fp1 = fopen("result.out", "w");
	double **a, *x, *t,*y;
	double X = PI, Y = PI, T = 0.2 ;
	a = (double**)malloc(NX * sizeof(double*));
	for (i = 0; i < NX; i++)  // цикл по строкам
		a[i] = (double*)malloc(NY * sizeof(double));
	x = (double*)malloc(NN * sizeof(double));
	y = (double*)malloc(NN * sizeof(double));
	t = (double*)malloc(NN * sizeof(double));
	//заполняем сетку
	hx = X / (NX-1);
	hy = Y / (NY-1);
	tau = T / (NT-1);
	for (int j = 0; j < NT; j++)
	{
		t[j] = tau* j ;
	}

	for (int k = 0; k < NX; k++)
	{
		x[k] = hx* k;
	}

	for (int k = 0; k < NT; k++)
	{
		y[k] = hy*k;
	}
	// гр.условия 


	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NY; j++)
		//a[i][j] = cos(3 * x[i]) * sin(4 * y[j]);		
			a[i][j] =exp((-pow((x[i]- PI / 2),2.0) - pow((y[j] - PI / 2), 2.0))/0.08);
	}
	
	
   for (int t = 1; t < NT; t++)
	 {
	   a = progonX(t, a);

	   a = progonY(t, a);
	   if (t == 35)
	   {
		   for (int i = 0; i < NX; i++)
		   {
			   for (int j = 0; j < NY; j++)

			   {
				   fprintf(fp1, "%lf %lf %lf\n", x[i], y[j], a[i][j]);
			   }
		   }
	   }
     }
   fclose(fp1);
}