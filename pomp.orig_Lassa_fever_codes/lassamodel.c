// -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>


//Ref: https://www.webpages.uidaho.edu/~barannyk/Teaching/LU_factorization_tridiagonal.pdf
//using LU decomposition of a tridiagonal matrix

void cubic_spline1(float xa[], float ya[], int n, float u[]);
void cubic_spline2(float xa[], float ya[], float u[], int n, float x, float *y);

// prototypes
double expit (double x);
double logit (double x);
void euler_multinomial (int ntrans, double n, double *rate, double *trans, double dt);
void negbin_rmeasure (int *n, int *index, double *X, double *cases);
void negbin_dmeasure (int *n, int *index, double *X, double *y1, double *y2, double *y3, double *f); // matching number of data
void basic_seir_b (double *Sh, double *Eh, double *Ih, double *Eq, double *Iq, double *H, double *R, double *C, 
                   double beta1, double beta2, double beta3, double muv, double muh, double sigma1, 
		   double gamma1, double tau1, double alpha1, double delta1, double sigma2, double tau2,
		   double alpha2, double delta2, double tau3, double delta3, double xi, 
                   double n, double rho, double dt);
void basic_sirs_pois (double *X, double *t1, double *t2, int *nstep,
			double *t_start, double *t_end,
		      int *nvar, int *np,
		      int *stateindex, int *parindex, double *temp1, double *temp2, double *temp3, int *seas_dim);

double expp(double gamma, double n1, double dt)   //gamma to infectious period in days
{
return(365.0*dt*n1/(1.0-exp(-n1*gamma*dt)));
}
double logg(double gamma, double n1, double dt)  //infectious period in days to gamma
{
return(-log(1.0-(365.0*dt*n1)/gamma)/(n1*dt));
}

void basic_seir_b (double *Sh, double *Eh, double *Ih, double *Eq, double *Iq, double *H, double *R, double *C, 
                   double beta1, double beta2, double beta3, double muv, double muh, double sigma1, 
		   double gamma1, double tau1, double alpha1, double delta1, double sigma2, double tau2,
		   double alpha2, double delta2, double tau3, double delta3, double xi, 
                   double n, double rho, double dt)
{
  double rate[20], trans[20], lambdah;


  if (!(R_FINITE(Sh[0]))) return;
  if (!(R_FINITE(Eh[0]))) return;
  if (!(R_FINITE(Ih[0]))) return;
  if (!(R_FINITE(Eq[0]))) return;
  if (!(R_FINITE(Iq[0]))) return;
  if (!(R_FINITE(H[0]))) return;
  if (!(R_FINITE(R[0]))) return;
  if (!(R_FINITE(C[0]))) return;

  if (!(R_FINITE(beta1))) return;
  if (!(R_FINITE(beta2))) return;
  if (!(R_FINITE(beta3))) return;
  if (!(R_FINITE(muv))) return;
  if (!(R_FINITE(muh))) return;
  if (!(R_FINITE(sigma1))) return;
  if (!(R_FINITE(gamma1))) return;
  if (!(R_FINITE(tau1))) return;
  if (!(R_FINITE(alpha1))) return;
  if (!(R_FINITE(delta1))) return;
  if (!(R_FINITE(sigma2))) return;
  if (!(R_FINITE(tau2))) return;
  if (!(R_FINITE(alpha2))) return;
  if (!(R_FINITE(delta2))) return;
  if (!(R_FINITE(tau3))) return;
  if (!(R_FINITE(delta3))) return;
  if (!(R_FINITE(xi))) return;
  if (!(R_FINITE(rho))) return;

// if rounding is not done elsewhere, we must do rounding here!
// the state variables are assumed to be positive integers!
//   S[0] = rint(S[0]);


  lambdah=(beta1*Ih[0]+beta2*Iq[0]+beta3*H[0]+muv)/n;
  rate[0] = lambdah; rate[1] = muh; rate[2] = sigma1; rate[3] = gamma1; rate[4] = muh;
  rate[5] = tau1; rate[6] = alpha1; rate[7] = delta1; rate[8] = muh; 
  rate[9] = sigma2; rate[10] = muh; rate[11] = tau2; rate[12] = alpha2; rate[13] = delta2;
  rate[14] = muh; rate[15] = tau3; rate[16] = delta3; rate[17] = muh; rate[18] = xi; rate[19] = muh;
  
  euler_multinomial(2,Sh[0],&rate[0],&trans[0],dt);
  euler_multinomial(3,Eh[0],&rate[2],&trans[2],dt);
  euler_multinomial(4,Ih[0],&rate[5],&trans[5],dt);
  euler_multinomial(2,Eq[0],&rate[9],&trans[9],dt);
  euler_multinomial(4,Iq[0],&rate[11],&trans[11],dt);
  euler_multinomial(3,H[0],&rate[15],&trans[15],dt);
  euler_multinomial(2,R[0],&rate[18],&trans[18],dt);
  Sh[0]  +=  (2500+trans[18] - trans[0] - trans[1]);
  Eh[0]  +=  (trans[0] - trans[2] - trans[3] - trans[4]);
  Ih[0]  +=  (trans[2] - trans[5] - trans[6] - trans[7] - trans[8]);
  Eq[0]  +=  (trans[3] - trans[9] - trans[10]);
  Iq[0]  +=  (trans[9] - trans[11]- trans[12]- trans[13]- trans[14]);
  H[0]   +=  (trans[6] + trans[12] - trans[15]- trans[16]-trans[17]);
  R[0]   +=  (trans[5] + trans[11] + trans[15]- trans[18]-trans[19]);
  C[0]   +=  (rint(rho*(trans[5] +trans[11])) + trans[6]+trans[12]+trans[7]+trans[13]+trans[16]);
//  if(Ih[0]<1)Ih[0]=1;
  if(Eh[0]+Ih[0]+Eq[0]+Iq[0]<1){Eh[0]+=1;Sh[0]-=1;}

}


#define BETA1      (x[parindex[0]])
#define BETA2      (x[parindex[1]])
#define BETA3      (x[parindex[2]])
#define MUH        (x[parindex[3]])
#define SIGMA1     (x[parindex[4]])
#define GAMMA1     (x[parindex[5]])
#define TAU1       (x[parindex[6]])
#define ALPHA1     (x[parindex[7]])
#define DELTA1     (x[parindex[8]])
#define SIGMA2     (x[parindex[9]])
#define TAU2       (x[parindex[10]])
#define ALPHA2     (x[parindex[11]])
#define DELTA2     (x[parindex[12]])
#define TAU3       (x[parindex[13]])
#define DELTA3     (x[parindex[14]])
#define XI         (x[parindex[15]])
#define RHO        (x[parindex[16]])
#define NM         (x[parindex[17]])
#define LOGMUV     (x[parindex[18]])
#define BMUV       (x[parindex[19]])
#define CMUV       (x[parindex[20]])
#define POP        (x[parindex[21]])
#define T0         (x[parindex[22]])

#define BSH      (x[stateindex[0]])
#define BEH      (x[stateindex[1]])
#define BIH      (x[stateindex[2]])
#define BEQ      (x[stateindex[3]])
#define BIQ      (x[stateindex[4]])
#define BH       (x[stateindex[5]])
#define BR       (x[stateindex[6]])
#define BCA      (x[stateindex[7]])
#define CSH      (x[stateindex[8]])
#define CEH      (x[stateindex[9]])
#define CIH      (x[stateindex[10]])
#define CEQ      (x[stateindex[11]])
#define CIQ      (x[stateindex[12]])
#define CH       (x[stateindex[13]])
#define CR       (x[stateindex[14]])
#define CCA      (x[stateindex[15]])
#define DSH      (x[stateindex[16]])
#define DEH      (x[stateindex[17]])
#define DIH      (x[stateindex[18]])
#define DEQ      (x[stateindex[19]])
#define DIQ      (x[stateindex[20]])
#define DH       (x[stateindex[21]])
#define DR       (x[stateindex[22]])
#define DCA      (x[stateindex[23]])


void basic_sirs_pois (double *X, double *t1, double *t2, int *nstep,
		      double *t_start, double *t_end,
		      int *nvar, int *np,
		      int *stateindex, int *parindex, double *temp1, double *temp2, double *temp3, int *seas_dim)
{
  double dt = (*t2 - *t1) / ((double) *nstep);
  float muv;
  double t, *x;
  float xa[*seas_dim] ,ya[*seas_dim], y2a[*seas_dim];
  int j, p, step;
  
  GetRNGstate();		// initialize R's pseudorandom number generator


  for (p = 0; p < *np; p++) {	// set the number of deaths to zero
    x = &X[*nvar * p];
    BCA = 0.0;
    CCA = 0.0;
    DCA = 0.0;
  }

   for (p = 0; p < *np; p++) {

   x = &X[*nvar * p];

	for(j = 0; j < *seas_dim; j++) {
	xa[j] = 0.0 + (1.0)/((float)*seas_dim-1.0)*j ;
	ya[j] = (float)(&LOGMUV)[j];
	}

  cubic_spline1(xa, ya, *seas_dim, y2a);
	if(p<0)
	for(j=0;j< *seas_dim;j++)
	Rprintf("%10.5f %10.5f %10.5f\n",xa[j+1],ya[j+1],y2a[j+1]);

  for (step = 0, t = *t1; step < *nstep; step++, t += dt) {
//  day=(int)floor(365.25*(t-2016.324-expit(T0)));
  cubic_spline2(xa, ya, y2a, *seas_dim, (t-*t_start)/(*t_end-*t_start), &muv);
//	muv=exp((&LOGMUV)[0]);
//	muv2=exp((&LOGMUV)[1]);
	if(p<0)
	Rprintf("%10.5f %10.5f\n",t,beta);

basic_seir_b(
   &BSH,  &BEH, &BIH, &BEQ, &BIQ, &BH, &BR, &BCA,
   exp(CMUV), exp(CMUV), exp(CMUV), exp(muv), exp(MUH), exp(SIGMA1), exp(GAMMA1), exp(TAU1), exp(ALPHA1), exp(DELTA1), exp(SIGMA2), exp(TAU2), 
   exp(ALPHA2), exp(DELTA2), exp(TAU3),exp(DELTA3), exp(XI), POP,  expit(RHO), dt);
basic_seir_b(
   &CSH,  &CEH, &CIH, &CEQ, &CIQ, &CH, &CR, &CCA,
   exp(CMUV), exp(CMUV), exp(CMUV), exp(muv), exp(MUH), exp(SIGMA1), exp(GAMMA1), exp(TAU1), exp(ALPHA1), exp(DELTA1), exp(SIGMA2), exp(TAU2), 
   exp(ALPHA2), exp(DELTA2), exp(TAU3),exp(DELTA3), exp(XI), POP,  expit(RHO), dt);
basic_seir_b(
   &DSH,  &DEH, &DIH, &DEQ, &DIQ, &DH, &DR, &DCA,
   exp(CMUV), exp(CMUV), exp(CMUV), exp(muv), exp(MUH), exp(SIGMA1), exp(GAMMA1), exp(TAU1), exp(ALPHA1), exp(DELTA1), exp(SIGMA2), exp(TAU2), 
   exp(ALPHA2), exp(DELTA2), exp(TAU3),exp(DELTA3), exp(XI), POP,  expit(RHO), dt);

    }
  }

  PutRNGstate();		// finished with R's random number generator

}


// simulate Euler-multinomial transitions
void euler_multinomial (int ntrans, double n, double *rate, double *trans, double dt) {
  double p = 0.0;
  int k;
  for (k = 0; k < ntrans; k++) p += rate[k]; // total event rate
  n = rbinom(n,1-exp(-p*dt));	// total number of events
  ntrans -= 1;
  for (k = 0; k < ntrans; k++) {
    trans[k] = ((n > 0) && (p > 0)) ? rbinom(n,rate[k]/p) : 0;
    n -= trans[k];
    p -= rate[k];
  }
  trans[ntrans] = n;
}


double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

double logit (double x) {
  return log(x/(1-x));
}

#define BC   (x[index[0]])
#define CC   (x[index[1]])
#define DC   (x[index[2]])
#define TAU  (x[index[3]])

void negbin_dmeasure (int *n, int *index, double *X, double *y1, double *y2,double *y3,double *f) {
  int p, nv = n[0], np = n[1];
  double *x,  size, prob, f1, f2, f3, tau, tol = 1.0e-18;
  for (p = 0; p < np; p++) {
    x = &X[nv*p];
    tau = exp(TAU);
    size = 1.0/tau;
    prob = 1.0/(1.0+BC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
     if (ISNA(*y1))
        f1 = 1.0;
        else
      f1 = dnbinom(*y1,size,prob,0)+tol;
    } else {
      f1 = tol;
    }
   prob = 1.0/(1.0+CC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
     if (ISNA(*y2))
        f2 = 1.0;
         else
      f2 = dnbinom(*y2,size,prob,0)+tol;
    } else {
      f2 = tol;
    }
   prob = 1.0/(1.0+DC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
     if (ISNA(*y3))
        f3 = 1.0;
         else
      f3 = dnbinom(*y3,size,prob,0)+tol;
    } else {
      f3 = tol;
    }
   	f[p]=f1*f2*f3;
  }
}

void negbin_rmeasure (int *n, int *index, double *X, double *cases) {
  int p, nv = n[0], np = n[1];
  double *x, tau, size, prob, tol = 1.0e-18;
  GetRNGstate();
  for (p = 0; p < np; p++) {
    x = &X[nv*p];
    tau = exp(TAU);
    size = 1.0/tau;
    prob = 1.0/(1.0+BC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
      cases[3*p] = rnbinom(size,prob+tol); // matching number of data
    } else {
      cases[3*p] = R_NaReal; // matching number of data
    }
	prob = 1.0/(1.0+CC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
      cases[3*p+1] = rnbinom(size,prob+tol);
    } else {
      cases[3*p+1] = R_NaReal;
    }
      prob = 1.0/(1.0+DC*tau);
    if (R_FINITE(size) && R_FINITE(prob)) {
      cases[3*p+2] = rnbinom(size,prob+tol);
    } else {
      cases[3*p+2] = R_NaReal;
    }
  }
  PutRNGstate();
}

#undef BC
#undef CC
#undef DC
#undef TAU



void cubic_spline1(float xa[], float ya[], int n, float u[])
{
    float A[n][n];
    float r[n];
    float h[n-1];
    float d[n-1];
    float yp1 = 0;
    float ypn = 0;
    float a[n], b[n], c[n];
    float b1[n], a1[n], r1[n];
    int i, j, k;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = 0;
        }
    }

    for (i = 0; i < n - 1; i++)
    {
        h[i] = xa[i + 1] - xa[i];
        d[i] = ya[i + 1] - ya[i];
    }

    for (i = 1; i < n - 1; i++)
    {
        A[i][i - 1] = h[i - 1];
        A[i][i] = 2 * (h[i - 1] + h[i]);
        A[i][i + 1] = h[i];
        r[i] = d[i] / h[i] - d[i - 1] / h[i - 1];
    }

    A[0][0] = 2 * h[0];
    A[0][1] = h[0];
    r[0] = d[0] / h[0] - yp1;


    A[n - 1][n - 2] = h[n - 2];
    A[n - 1][n - 1] = 2 * h[n - 2];
    r[n - 1] = ypn - d[n - 2] / h[n - 2];

    for (i = 0; i < n; i++)
    {
        r[i] *= 6;
    }

    for (j = 0; j < n; j++)
    {
        a[j] = 0;
        if (j > 0)
        {
            a[j] = A[j][j - 1];
        }
        b[j] = A[j][j];
        c[j] = 0;
        if (j < n - 1)
        {
            c[j] = A[j][j + 1];
        }
    }

    for (i = 0; i < n; i++)
    {
        b1[i] = b[i];
        a1[i] = a[i];
        r1[i] = r[i];
    }

    b1[0] = b[0];
    for (k = 1; k < n; k++)
    {
        a1[k] = a[k] / b1[k - 1];
        b1[k] = b[k] - a[k] / b1[k - 1] * c[k - 1];
    }

    r1[0] = r[0];
    for (k = 1; k < n; k++)
    {
        r1[k] = r[k] - a1[k] * r1[k - 1];
    }

    u[n - 1] = r1[n - 1] / b1[n - 1];
    for (k = n - 2; k >= 0; k--)
    {
        u[k] = (r1[k] - c[k] * u[k + 1]) / b1[k];
    }

}


void cubic_spline2(float xa[], float ya[], float u[], int n, float x, float *y){
    int o = 0, i = n-1, k;

    while (1) {
        if (i - o <= 1) break;
        k = round((i + o) / 2);
        if (xa[k] > x) i = k;
        else o = k;
    }

    *y = ((xa[i] - x) / (xa[i] - xa[o])) * ya[o] + ((x - xa[o]) / (xa[i] - xa[o])) * ya[i] +
         ((((xa[i] - x) / (xa[i] - xa[o])) * ((xa[i] - x) / (xa[i] - xa[o])) * ((xa[i] - x) / (xa[i] - xa[o])) - ((xa[i] - x) / (xa[i] - xa[o]))) * u[o] +
          (((x - xa[o]) / (xa[i] - xa[o])) * ((x - xa[o]) / (xa[i] - xa[o])) * ((x - xa[o]) / (xa[i] - xa[o]))-((x - xa[o]) / (xa[i] - xa[o]))) * u[i])*((xa[i] - xa[o])*(xa[i] - xa[o])) / (3.0 * 2);
}


