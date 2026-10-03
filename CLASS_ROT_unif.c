/************************************************************************** 
*                                                                         *
**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"

char version[30] = "1.0, September 2015";


 #define PI 3.1415926535              /* what else */
 #define SQ(x) ((x)*(x))              /* square macro */
 #define enne 1.5                     /* polytropic index */
 #define step 1.05                    /* step in axis ratio to build the sequence of models */
 #define MODNUM_MAX 100
 #define MAXITS 200
 #define NT 8                         /* grid angular dimension */
 #define NR 22                         /* grid radial dimension */
 #define DNMAX 11                     /* number of terms in Legendre polynomial expansions */
 #define TINY 1.0e-20                 /* A small number */
 #define EPS 1.0e-4
 #define TOLX 1.0e-4
 #define TOLF 1.0e-7
 #define ALF 1.0e-4                   /* Ensures sufficient decrease in function value. */
 #define STPMX 100.0
 #define TOLMIN 1.0e-6

 #define FREERETURN {free_dvector(fvec,1,n);free_dvector(xold,1,n);free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);free_ivector(indx,1,n);return;}

/* Here MAXITS is the maximum number of iterations; 
   TOLF sets the convergence criterion on function values; 
   TOLMIN sets the criterion for deciding whether spurious convergence to a minimum of FNORM has occurred; 
   TOLX is the convergence criterion on δx; 
   STPMX is the scaled maximum step length allowed in line searches. */

void load_sphere( char sphere[]);
void make_grid(double y[]);
void usrfun(double x[],int n,double *fvec,double **fjac);
void newt(double x[], int n, int *check,void (*vecfunc)(int, double [], double []));
void vecfunc(int n,double x[],double f[]);
void print_header(void);
void (*nrfuncv)(int n, double v[], double f[]);
double FNORM(double x[]);
double plgndr(int l, int m, double x);

char file_name[80] = "no sphere file specified", 
     eos_type[30] = "tab"; 

int ITERS,
    CKC=0,
    nn,
    a_check=0,                       /* indicates if iteration diverged */ 
    print_option,                    /* select print out */
    num_it,                          /* counter of number of iterations in N-R method */
    modnum=0,                        /* model number */
    cntr = 0;                        /* control parameter for the model: - 0 by default
                                                                         - 1 problems by Nwt-Raph routine 
                                                                         - 2 no Hoiland crituerium
                                                                         - 3 eq eff grav >0
                                                                         - 4 no virial theorem
                                                                         - 5 sequence exceeds MODNUM_MAX models 
                                                                         - 6 OK */


double obl;                      /* initial value for axis ratio defined as r_eq/r_pol */
double *SPH_val;         /* spherical model */
double *init;            /* variables for the guess model */
double *theta;                  /* grid points in angular direction */
double *THETA;                  /* coefficient for integration in angular direction */
double **ERRE;               /* coefficient for integration in radial direction */
double *LEG_POL_2c_i;     /* vector of Legendre Polynomials used for gravitational potential */
double *fgr;       /* vector of Green's functions coefficients */
double Fsys_x[NT*(NR+1)+3];          /* vector which contains value of system functions computed in a certain point x */
double a;                             /* parameter for rotation law in system = a*r_s, where r_s is the radius of sphere */
double aa[3];                         /* dimensionless parameter for differential rotation */                                                                        
double *fvec;                                                                         

/********************************************************************/
/* Routine for LU decomposition from Numerical Recipies in C        */
/********************************************************************/

void ludcmp(double **a, int n, int *indx, double *d)

/* Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
   permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
   indx[1..n] is an output vector that records the row permutation effected by the partial
   pivoting; d is output as ±1 depending on whether the number of row interchanges was even
   or odd, respectively. This routine is used in combination with lubksb to solve linear equations
   or invert a matrix. */
{

int i,
    imax,
    j,
    k;

double big,
      dum,
      sum,
      temp;

double *vv;                                        // vv stores the implicit scaling of each row.

vv = dvector(1,n);
*d=1.0;                                           // No row interchanges yet.

for (i=1;i<=n;i++) {                              // Loop over rows to get the implicit scaling information.
    big=0.0;

    for (j=1;j<=n;j++) {
         if ((temp=fabs(a[i][j])) > big) big=temp;
       }
    
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");  // No nonzero largest element.
                                                                   
    vv[i]=1.0/big;                                                 // Save the scaling.
    
}

for (j=1;j<=n;j++) {                                               // This is the loop over columns of Crout’s method.
 
 for (i=1;i<j;i++) {                                               // This is equation (2.3.12) except for i = j.
      sum=a[i][j];
      
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      
      a[i][j]=sum;
  }

 big=0.0;                                                          // Initialize for the search for largest pivot element.

 for (i=j;i<=n;i++) {                                              // This is i = j of equation (2.3.12) and i = j + 1 ...N
    sum=a[i][j];                                                      // of equation (2.3.13).

  for (k=1;k<j;k++) sum -= a[i][k]*a[k][j];
  
  a[i][j]=sum;
 
  if ( (dum=vv[i]*fabs(sum)) >= big) {
                                                              // Is the figure of merit for the pivot better than the best so far?
    big=dum;
    imax=i;
  }
 
 }

 if (j != imax) {                                              // Do we need to interchange rows?

  for (k=1;k<=n;k++) {                                          // Yes, do so...
     dum=a[imax][k];
     a[imax][k]=a[j][k];
    a[j][k]=dum;
  }
 
  *d = -(*d);                                                   // ...and change the parity of d.
  vv[imax]=vv[j];                                               // Also interchange the scale factor.
 }

 indx[j]=imax;
 
 if (a[j][j] == 0.0) a[j][j]=TINY;
                                     // If the pivot element is zero the matrix is singular (at least to the precision of the
                                     // algorithm). For some applications on singular matrices, it is desirable to substitute
                                     // TINY for zero.

 if (j != n) {                       // Now, finally, divide by the pivot element.

   dum=1.0/(a[j][j]);
   
   for (i=j+1;i<=n;i++) a[i][j] *= dum;
  }

}                                     // Go back for the next column in the reduction.

free_dvector(vv,1,n);

return;

}

/********************************************************************/
/* Routine to solve linear system from Numerical Recipies in C      */
/********************************************************************/

void lubksb(double **a, int n, int *indx, double b[])

/* Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
   A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
   as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
   B, and returns with the solution vector X. a, n, and indx are not modified by this routine
   and can be left in place for successive calls with different right-hand sides b. This routine takes
   into account the possibility that b will begin with many zero elements, so it is efficient for use
   in matrix inversion. */
{

 int i,
     ii=0,
     ip,
     j;

 double sum;

 for (i=1;i<=n;i++) {                             // When ii is set to a positive value, it will become the
                                                  // index of the first nonvanishing element of b. We now
                                                  // do the forward substitution, equation (2.3.6). The
                                                  // only new wrinkle is to unscramble the permutation
                                                  // as we go.
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    
    if (ii)  
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;                           // A nonzero element was encountered, so from now on we
                                                  // will have to do the sums in the loop above.
    b[i]=sum;                                      
  }
 
 for (i=n;i>=1;i--) {                              // Now we do the backsubstitution, equation (2.3.7).
     sum=b[i];
     
     for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
     b[i]=sum/a[i][i];                             // Store a component of the solution vector X.
  
  }                                                // All done!
 
 return;
}

/***********************************************************************************************/
/* Routine for linear search used in Newton-Raphson method from Numerical Recipies in C        */
/***********************************************************************************************/

void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[], double *f, double stpmax, int *check, double (*func)(double []))

/* Given an n-dimensional point xold[1..n], the value of the function and gradient there, 
   fold and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction 
   p from xold where the function func has decreased “sufficiently.” The new function value 
   is returned in f. stpmax is an input quantity that limits the length of the steps so that 
   you do not try to evaluate the function in regions where it is undefined or subject to overflow. 
   p is usually the Newton direction. The output quantity check is false (0) on a normal exit. 
   It is true (1) when x is too close to xold. In a minimization algorithm, this usually signals 
   convergence and can be ignored. However, in a zero-finding algorithm the calling program 
   should check whether the convergence is spurious. Some “difficult” problems may require double 
   precision in this routine. */

   
{

int i;

double a,
      alam,
      alam2,
      alamin,
      b,
      disc,
      f2,
      rhs1,
      rhs2,
      slope,
      sum, 
      temp,
      test,
      tmplam;

*check=0;

for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i]; 

sum=sqrt(sum);

if (sum > stpmax) 
  for (i=1;i<=n;i++) p[i] *= stpmax/sum;                    // Scale if attempted step is too big.

for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];

if (slope >= 0.0) nrerror("Roundoff problem in lnsrch."); 

test=0.0;                                                   // Compute λmin. 

for (i=1;i<=n;i++) {
  temp=fabs(p[i])/DMAX(fabs(xold[i]),1.0);
  if (temp > test) test=temp; 
}
alamin=TOLX/test;
alam=1.0;                                                   // Always try full Newton step first.

for(;;) {                                                   // Start of iteration loop.
   for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i]; 
   *f=(*func)(x);
   if (alam < alamin) {
   for (i=1;i<=n;i++) x[i]=xold[i]; 
   *check=1;
   return;
   }
   else if (*f <= fold+ALF*alam*slope) return; 
   else {
       if (alam == 1.0)
          tmplam = -slope/(2.0*(*f-fold-slope));
       else {
          rhs1 = *f-fold-alam*slope;
          rhs2=f2-fold-alam2*slope; 
          a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2); 
          b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2); 
          if (a == 0.0) tmplam = -slope/(2.0*b);
          else {
              disc=b*b-3.0*a*slope;
              if (disc < 0.0) tmplam=0.5*alam;
              else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a); 
              else tmplam=-slope/(b+sqrt(disc));
               }
          if (tmplam > 0.5*alam)
              tmplam=0.5*alam;
            } 
        }
   alam2=alam;
   f2 = *f; 
   alam=DMAX(tmplam,0.1*alam);                               // λ ≤ 0.5λ1.

 }
 
 return;
}

/*****************************************************************************************/
/* Globally convergent Newton-Raphson routine from Numerical Recipies in C               */
/*****************************************************************************************/

void newt(double x[], int n, int *check,void (*vecfunc)(int, double [], double []))

/* Given an initial guess x[1..n] for a root in n dimensions, find the root by 
   a globally convergent Newton’s method. The vector of functions to be zeroed, 
   called fvec[1..n] in the routine below, is returned by the user-supplied routine 
   vecfunc(n,x,fvec). The output quantity check is false (0) on a normal return and 
   true (1) if the routine has converged to a local minimum of the function FNORM 
   defined below. In this case try restarting from a different initial guess.*/

{

void fdjac(int n, double x[], double fvec[], double **df,
void (*vecfunc)(int, double [], double []));
double FNORM(double x[]);
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],double *f, double stpmax, int *check, double (*func)(double [])); 
void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);

int i,
    its,
    j,
    *indx;

double d,
      den,
      f,
      fold,
      stpmax,
      sum,
      temp,
      test,
      **fjac,
      *g,
      *p,
      *xold;

indx=ivector(1,n); 
fjac=dmatrix(1,n,1,n); 
g=dvector(1,n); 
p=dvector(1,n); 
xold=dvector(1,n); 
fvec=dvector(1,n);                                   // Define global variables.
nn=n; 
nrfuncv=vecfunc; 
f=FNORM(x);                                          // fvec is also computed by this call.

test=0.0;                                           // Test for initial guess being a root. 
                                                    // Use more stringent test than simply TOLF
for (i=1;i<=n;i++)
  if (fabs(fvec[i]) > test) test=fabs(fvec[i]);

if (test < 0.01*TOLF) { 
  *check=0;
  FREERETURN 
}

for (sum=0.0,i=1;i<=n;i++) sum += SQ(x[i]);        // Calculate stpmax for line searches.

stpmax=STPMX*DMAX(sqrt(sum),(double)n);

for (its=1;its<=MAXITS;its++) {                     // Start of iteration loop.

   ITERS = its;

   fdjac(n,x,fvec,fjac,vecfunc);
 
   // If analytic Jacobian is available, you can replace the routine fdjac below with your own routine.
   
   for (i=1;i<=n;i++) {                             // Compute ∇f for the line search.
       for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
       g[i]=sum; 
   }

   for (i=1;i<=n;i++) xold[i]=x[i];                // Store x
   fold=f;                                         // and f
   for (i=1;i<=n;i++) p[i] = -fvec[i];             // Right hand side for linear equation
   
   ludcmp(fjac,n,indx,&d);                         // Solve linear equation by LU decomposition
   lubksb(fjac,n,indx,p);
   lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,FNORM); 
                    // lnsrch returns new x and f. It also calculates fvec at the new x when it calls FNORM.

   test=0.0;                                       // Test for convergence on function values.

   for (i=1;i<=n;i++) 
       if (fabs(fvec[i]) > test) test=fabs(fvec[i]); 

   if (test < TOLF) {
   *check=0;
   FREERETURN 
   }

   if (*check) {                                  // Check for gradient of f zero, i.e., spurious convergence.
       test=0.0;
       den=DMAX(f,0.5*n); 
       for (i=1;i<=n;i++) {
           temp=fabs(g[i])*DMAX(fabs(x[i]),1.0)/den;
           if (temp > test) test=temp; 
         }
       *check=(test < TOLMIN ? 1 : 0);
       FREERETURN 
     }

   test=0.0;                                     // Test for convergence on δx.
   for (i=1;i<=n;i++) {
       temp=(fabs(x[i]-xold[i]))/DMAX(fabs(x[i]),1.0);
       if (temp > test) test=temp; 
     }

   if (test < TOLX) FREERETURN 
 }
   nrerror("MAXITS exceeded in newt"); 

//return its;
}

/***********************************************************************************************/
/* Routine to compute Jacobian of a vector function from Numerical Recipies in C               */
/***********************************************************************************************/

void fdjac(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double []))

/* Computes forward-difference approximation to Jacobian. On input, x[1..n] is the point at
   which the Jacobian is to be evaluated, fvec[1..n] is the vector of function values at the
   point, and vecfunc(n,x,f) is a user-supplied routine that returns the vector of functions at
   x. On output, df[1..n][1..n] is the Jacobian array. */
 {

 int i,
     j;

 double h,
       temp,
       *f;
 
 f = dvector(1,n);
 
 for (j=1;j<=n;j++) {

     temp=x[j];
     h=EPS*fabs(temp);
    
     if (h == 0.0) h=EPS;
     
     x[j]=temp+h;                          // Trick to reduce finite precision error.
     h=x[j]-temp;
     
     (*vecfunc)(n,x,f);
     
     x[j]=temp;
     for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;    // Forward difference formula.

  } 
  
  free_dvector(f,1,n);
  return;
}

/************************************************************************************/
/* Routine for norm of a vector function from Numerical Recipies in C               */
/************************************************************************************/

//extern int nn;
//extern double *fvec;
//extern void (*nrfuncv)(int n, double v[], double f[]);

double FNORM(double x[])

/* Returns f = 1/2 F · F at x. The global pointer *nrfuncv points to a routine that returns the
vector of functions at x. It is set to point to a user-supplied routine in the calling program. 
Global variables also communicate the function values back to the calling program. */

{
  int i;
  double sum;

  (*nrfuncv)(nn,x,fvec);
  for (sum=0.0,i=1;i<=nn;i++) sum += SQ(fvec[i]); 

  return 0.5*sum;
}
       
/*********************************************************************/
/* Create computational mesh for weight for integration              */
/* Points in the angular direction are stored in the array theta[i]. */
/* Points in the r-direction are stored in the array s_pg[j].        */
/*********************************************************************/
void make_grid(double y[])                                    
{ 
  int i,          /* counter theta */
      j;          /* counter r */

 for(i=1;i<=NT;i++) {
     if(i==0){
      THETA[i] = (PI/(4.0*(NT-1.0)));
      }
     else if(i==NT) {
      THETA[i] = (PI/(4.0*(NT-1.0)));
     }
     else {
      THETA[i] = (PI/(2.0*(NT-1.0)));
     }
 }
 
 /* Grid points in r directions are x[i]*j/NR, for i=1..NT and j=1..NR, x[i] = r_surf[i] */

 for(i=1;i<=NT;i++) {
  for(j=1;j<=NR;j++) {
    if(j==1) {
      ERRE[i][j] = (y[i])/(3.0*NR);
    }
    else if(j==(NR)) {
      ERRE[i][j] = (y[i])/(3.0*NR);         
    }                                    
    else if((j % 2)==0) {
      ERRE[i][j] = (4.0*y[i])/(3.0*NR);  
    }                                   
    else {
      ERRE[i][j] = (2.0*y[i])/(3.0*NR);
    }
  }
 }
 
 return;
}

/***********************************************************************************************/
/* Routine to Legendre Polynomials of given theta. The first is taken from Numerical Recipies  */
/***********************************************************************************************/

double plgndr(int l, int m, double x)
/* Computes the associated Legendre polynomial P^l_m(x). Here m and l are integers
   satisfying 0 ≤ m ≤ l, while x lies in the range −1 ≤ x ≤ 1. */
{

void nrerror(char error_text[]); 

double fact,
       pll,
       pmm,
       pmmp1,
       somx2; 

int i,
    ll;

if (m < 0 || m > l || fabs(x) > 1.0) nrerror("Bad arguments in routine plgndr.\n");
pmm=1.0; 

if (m > 0) {                                   //Compute P_m^m. 
   somx2=sqrt((1.0-x)*(1.0+x)); 
   fact=1.0;
   for (i=1;i<=m;i++) {
    pmm *= -fact*somx2;
    fact += 2.0; 
    }
}

if (l == m)
    return pmm;

else {                                         //Compute P_m^m+1
    pmmp1=x*(2*m+1)*pmm; 
    
    if (l == (m+1))
       return pmmp1;
    
    else {                                     //Compute P_l^m, l > m+1.
       for (ll=m+2;ll<=l;ll++) {
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m); pmm=pmmp1;
           pmmp1=pll;
           }
       return pll; 
     }
    }   

}


void leg_polyn(){

/* Order is LEG_POL_2c_i[i + c*NT] = P(2*c,cos(theta[i])) */
     
     double x;

     int c,
         i;

     for (i=1;i<=NT;i++){
      x=cos(theta[i]);
      for (c=0;c<DNMAX;c++){
        LEG_POL_2c_i [ i + c*NT ] = plgndr((2*c),0,x);
      }
      
     }
     return;
}

/*************************************************************************/
/* Load SPHERE file. Format:  NT*(NR+1)+10 vector                         */ 
/*************************************************************************/

void load_sphere(char sphere[])
{
 int i,                         /* counter */
     n_var;                     /* number of variables */

 double y;         /* spherical values ordered as listed below */
        
 FILE *f_sphere;                /* pointer to sphere */
 

    /* OPEN FILE TO READ */

    if((f_sphere=fopen(sphere,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",sphere); 
       exit(0);
    }

 
    /* COMPUTE NUMBER OF VARIABLES AND PHYSICAL QUANTITIES */

    n_var = NT*(NR+1)+11;


    /* READ SPHERE INPUT */

    for(i=1;i<n_var;i++) {  
       
/*************************************** HERE **************************************************************************/

       fscanf(f_sphere,"%lf ",&y) ;	

/*  the order is the following: - NT surface radii
                                -	NT*NR values of dimensionless density following radii gridpoints from center to surface, from pole to equator
                                - 1 dimensionless central angular velocity ( \nu_c = \Omega^2/(4*\Pi*G*\rho_c) )
                                - 1 dimensionless gravitational potential at the pole
                                -------------------------------------------------------------- NT*(NR+1)+2
                                - dimensionless kinetic energy
                                - dimensionless gravitational potential energy
                                - dimensionless thermal energy
                                - dimensionless total mass
                                - dimensionless total angular momentum
                                - dimensionless surface effective gravity at the equator
                                - cntr value to check if the sequence should be stopped
                                - VT = | ( 2*T + W + (3/n)*U ) / W |
                                --------------------------------------------------------------- NT*(NR+1)+10 */
        SPH_val[i] = y;
       }

    return; 
}

/***********************************************************/
/* Routine for rotation law in cylindrical coordinates     */
/***********************************************************/

// NB: xi_cil = r * sin(theta), thus OMEGA(xi_cil) = OMEGA(r,theta) if one considers spherical coordinates

// The routine returns the value of angular velocity at a certain distance from the axis 

double DIFF_ROT_LAW(double x[])
{
  double OM;

  OM = sqrt(x[NT*(NR+1)+1]);
  
  return OM;
}

/*****************************************/
/* Routine to compute Green's Functions  */
/*****************************************/

void FGREEN(double x[]){
  
  int i,
      j,
      p,
      q,
      c,
      v;

  for (i=1;i<=NT;i++){
    for (j=1;j<=NR;j++){
      for (p=1;p<=NT;p++){
        for (q=1;q<=NR;q++){
          for (c=0;c<DNMAX;c++){
            v = i + (j-1)*NT + (p-1)*NT*NR + (q-1)*NT*NT*NR + c*NR*NR*NT*NT;
            if ((x[p]*q/NR) > (x[i]*j/NR)) fgr[ v ] = pow((x[i]*j/NR),(2.0*c))/pow((x[p]*q/NR),(2.0*c+1.0));
            else fgr[ v ] = pow((x[p]*q/NR),(2.0*c))/pow((x[i]*j/NR),(2.0*c+1.0));
          }
        }
      }
    }
  }

  return;
}

/******************************************************************************************************************************/
/* Routine that returns the values of functions in the system of equations f[1,..,n], calculated at a certain point x[1,..,n] */
/******************************************************************************************************************************/

void vecfunc(int n,double x[],double f[]){

void FGREEN(double x[]);

int i,
    j,
    v,
    p,
    q,
    c; 

double GPC,
       GPG[NT+1][NR+1],
       CPC,
       CPG[NT+1][NR+1];

double sum_pot_r = 0.0,     // intermediate sum in r direction
       sum_pot_th = 0.0,    // intermediate sum in theta direction
       sum_pot_c = 0.0;     // intermediate sum in n

double temp=0.0,
       temp_in=0.0;

FGREEN(x);

  for (p=1;p<=NT;p++){
    temp_in = 0.0;
    for (q=1;q<=NR;q++){
      temp_in = temp_in + ( ((x[p]*q)/NR)* ERRE[p][q] * pow(x[NT+(p-1)*NR+q],enne) );
    }
    temp = temp + (sin(theta[p]) * THETA[p] * temp_in);
  }

  GPC = - temp;

  for (i=1;i<=NT;i++)
  {
    for (j=1;j<=NR;j++)
    {
      sum_pot_th = 0.0;          
      for(p=1;p<=NT;p++)
      {
        sum_pot_r = 0.0;
        for (q=1;q<=NR;q++)
        {
          sum_pot_c = 0.0; 
          for (c=0;c<DNMAX;c++)
          {
             sum_pot_c +=  fgr[ i + (j-1)*NT + (p-1)*NT*NR + (q-1)*NT*NT*NR + c*NR*NR*NT*NT ] * LEG_POL_2c_i[ i + c*NT ] * LEG_POL_2c_i[ p + c*NT ] * pow(x[ NT+((p-1)*NR)+q ],enne);
          }
          sum_pot_r +=  SQ((x[p]*q/NR))*ERRE[p][q]*sum_pot_c;
        }
        sum_pot_th += sin(theta[p])*THETA[p]*sum_pot_r;
      }
      GPG[i][j] = - sum_pot_th; 
    }
  }

/*CPC = (SQ(a)*x[NT*(NR+1)+1])/2.0;*/

for (i=1;i<=NT;i++) {
  for (j=1;j<=NR;j++){
    CPG[i][j] = - 0.5*(x[NT*(NR+1)+1]*SQ((x[i]*sin(theta[i])*j/NR)));
  }
}


f[1] = 1.0 + GPC /*+ CPC*/ - x[NT*(NR+1)+2];

v = 1;

for (i=1;i<=NT;i++){
  for(j=1;j<=NR;j++){
    v +=1;
    f[v] = x[NT+NR*(i-1)+j] + GPG[i][j] + CPG[i][j] - x[NT*(NR+1)+2];  
  }
}

for (i=1;i<=NT;i++) f[1+NT*NR+i] = x[NT+(i-1)*NR+NR];

if(n==(NT*(NR+1)+2)) f[n] = x[NT]-(obl*x[1]);
else printf("ERROR in SYSTEM DIMENSION's definition\n");

return;
}

/***********************************************************************************/
/* Routine to compute physical quantities stored in x[NT*(NR+1)+3,..,NT*(NR+1)+10] */
/***********************************************************************************/

/* x[NT*(NR+1)+3] - dimensionless kinetic energy
   x[NT*(NR+1)+4] - dimensionless gravitational potential energy
   x[NT*(NR+1)+5] - dimensionless thermal energy, defined as n/(n+1) * int(sin(theta)*dtheta*int(sigma^(n+1)*xi^2,xi=0..xi_surf(theta)),theta=0..PI/2), here being xi the dimensionless spherical radius
   x[NT*(NR+1)+6] - dimensionless total mass
   x[NT*(NR+1)+7] - dimensionless total angular momentum
   x[NT*(NR+1)+8] - dimensionless surface effective gravity at the equator 
   x[NT*(NR+1)+9] - cntr value for the model
   x[NT*(NR+1)+10]- V.T. */

void PHYS_QUANT(double x[], int modn){
  int p,
      q;
  double temp_p = 0.0,
        temp_q = 0.0;
  void CONTROL_MOD(double x[],int modn);
  
  for(p=1;p<=NT;p++){
    temp_q = 0.0;
    for(q=1;q<=NR;q++){
      temp_q += pow((x[p]*q/NR),4.0) * ERRE[p][q] * pow((x[NT+(p-1)*NR+q]),enne) * SQ(DIFF_ROT_LAW(x));
    }
    temp_p += pow(sin(theta[p]),3.0)* THETA[p] * temp_q;
  }

  x[NT*(NR+1)+3] = 0.5 * temp_p;

  temp_p = 0.0;
  temp_q = 0.0;

  for(p=1;p<=NT;p++){
    temp_q = 0.0;
    for(q=1;q<=NR;q++){
      temp_q += SQ((x[p]*q/NR)) * ERRE[p][q] * pow((x[NT+(p-1)*NR+q]),enne) * ( x[NT*(NR+1)+2] - x[NT+(p-1)*NR+q] + 0.5*x[NT*(NR+1)+1]*SQ(((x[p]*q/NR)*sin(theta[p]))) );
    }
    temp_p += sin(theta[p]) * THETA[p] * temp_q;
  }

  x[NT*(NR+1)+4] = 0.5 * temp_p;

  temp_p = 0.0;
  temp_q = 0.0;

  for(p=1;p<=NT;p++){
    temp_q = 0.0;
    for(q=1;q<=NR;q++){
      temp_q += SQ((x[p]*q/NR)) * ERRE[p][q] * pow((x[NT+(p-1)*NR+q]),enne+1.0);
    }
    temp_p += sin(theta[p]) * THETA[p] * temp_q;
  }

  x[NT*(NR+1)+5] = enne/(enne+1.0) * temp_p;

  temp_p = 0.0;
  temp_q = 0.0;

  for(p=1;p<=NT;p++){
    temp_q = 0.0;
    for(q=1;q<=NR;q++){
      temp_q += SQ((x[p]*q/NR)) * ERRE[p][q] * pow((x[NT+(p-1)*NR+q]),enne);
    }
    temp_p += sin(theta[p]) * THETA[p] * temp_q;
  }

  x[NT*(NR+1)+6] = temp_p;

  temp_p = 0.0;
  temp_q = 0.0;
  
  for(p=1;p<=NT;p++){
    temp_q = 0.0;
    for(q=1;q<=NR;q++){
      temp_q += pow((x[p]*q/NR),4.0) * ERRE[p][q] * pow((x[NT+(p-1)*NR+q]),enne) * DIFF_ROT_LAW(x);
    }
    temp_p += pow(sin(theta[p]),3.0)* THETA[p] * temp_q;
  }

  x[NT*(NR+1)+7] = temp_p;

  temp_p = 0.0;
  temp_q = 0.0;

  /* Using a three-points approximation for sigma(xi,theta) it turns out that
     ((d sigma)/(d xi))(xi_surf,Pi/2) = (1/2)*(NR/xi_surf)*( 3*sigma(xi_surf,Pi/2) +
     - 4*sigma(xi_surf*(NR-1)/NR,Pi/2) + sigma(xi_surf*(NR-2)/NR,Pi/2) ) */

  x[NT*(NR+1)+8] = 0.5 * (NR/x[NT]) * ( 3.0*x[NT+(NT-1)*NR+NR] - 4.0*x[NT+(NT-1)*NR+NR-1] + x[NT+(NT-1)*NR+NR-2] );

  CONTROL_MOD(x,modn);

  x[NT*(NR+1)+9] = cntr;
  
  x[NT*(NR+1)+10] = fabs(( (2.0*x[NT*(NR+1)+3]) + x[NT*(NR+1)+4] + ((3.0/enne)*x[NT*(NR+1)+5]) )/(x[NT*(NR+1)+4]));
  
  return;
}

/*******************************************/
/* Routine to assign cntr value for model  */
/*******************************************/

void CONTROL_MOD(double x[],int modn){
  
  int i,
      isnegative=0,
      temp;

  double part_diff_omega=0.0,
         part_diff_r=0.0,
         DIFF_ANG_MOM=0.0;

  /*for (i=1;i<=NT;i++){
      part_diff_omega=0.0;
      part_diff_r=0.0;
      DIFF_ANG_MOM=0.0;
      if(i<=2){part_diff_omega = 0.5 * ((NT-1.0)/(PI/2.0)) * ( -3.0*(DIFF_ROT_LAW((x[i]*sin(theta[i])),x)*SQ(x[i]*sin(theta[i]))) + 4.0*(DIFF_ROT_LAW((x[i+1]*sin(theta[i+1])),x)*SQ(x[i+1]*sin(theta[i+1]))) - (DIFF_ROT_LAW((x[i+2]*sin(theta[i+2])),x)*SQ(x[i+2]*sin(theta[i+2]))) );}
      else{part_diff_omega = 0.5 * ((NT-1.0)/(PI/2.0)) * ( 3.0*(DIFF_ROT_LAW((x[i]*sin(theta[i])),x)*SQ(x[i]*sin(theta[i]))) - 4.0*(DIFF_ROT_LAW((x[i-1]*sin(theta[i-1])),x)*SQ(x[i-1]*sin(theta[i-1]))) + (DIFF_ROT_LAW((x[i-2]*sin(theta[i-2])),x)*SQ(x[i-2]*sin(theta[i-2]))) );}
                            
      /* NB: in part_diff_omega the three points approx is used with
         forward substituion if i<=2 or backward substitution
         otherwise */

     /* part_diff_r = 0.5 * ((NR)/(x[i])) * ( 3.0*(DIFF_ROT_LAW((x[i]*sin(theta[i])),x)*SQ(x[i]*sin(theta[i]))) - 4.0*(DIFF_ROT_LAW((x[i]*((NR-1.0)/NR)*sin(theta[i])),x)*SQ(x[i]*((NR-1.0)/NR)*sin(theta[i]))) + (DIFF_ROT_LAW((x[i]*((NR-2.0)/NR)*sin(theta[i])),x)*SQ(x[i]*((NR-2.0)/NR)*sin(theta[i]))) );*/
                              
      /* NB: in partial differentiation with respect to r always
         three points approximation with backard substitution
         is adopted, as it is computed always in the surface
         radius */

     /* DIFF_ANG_MOM = part_diff_r * (sin(theta[i])) + part_diff_omega*(cos(theta[i])/x[i]);*/

      /* chain rule for multivariable calculus:
         df/dxi = df/dr * sin(theta) + df/dtheta * (cos(theta)/r),
         where (xi,z,phi) are cylindrical coordinates and
         (r,theta,phi) are spherical ones */
    
    /*  if(DIFF_ANG_MOM<0.0) temp=1;
      else temp=0;

      isnegative += temp;
      
    }*/


  if(ITERS>=MAXITS){cntr = 1;
                    return;}
  
  /*else if(isnegative > 0) {cntr = 2;
                           return;}*/
             
  else if((x[NT*(NR+1)+8])>0) {cntr=3;
                             return;}
  else if((x[NT*(NR+1)+3]/fabs(x[NT*(NR+1)+4]))-0.5>0.0){cntr=4;
                                               return;}
  else if(modn - MODNUM_MAX > 0) {cntr=5;
                                  return;}
  else if(modn - MODNUM_MAX <=0) {cntr=6;
                                  return;}
  
  
  return;
}


/*******************************************/
/* Routine to print header                 */
/*******************************************/

void print_header(void){

 printf("\n");
 printf("------------------------------------------------------------------------------------\n");   
 printf("NTxNR=%dx%d, a=%lf  accuracy in x=%1.0e, accuracy in f=%1.0e\n",NT,NR,a,TOLX,TOLF);

 printf("------------------------------------------------------------------------------------\n");   
 printf("   NT*r_surf[i]   |    (NT*NR)*rho[i,j]    |   nu_c   |   Phi_g_pol   |   E_kin   |  E_pot  |  U  |  M   |  J   |  eq effgrav  |  cntr  |  V.T.\n");

 printf("------------------------------------------------------------------------------------\n");   
 
 return;
}

/*******************************************/
/*      MAIN CYCLE                         */
/*******************************************/

int main(int argc, char *argv[])
{ 
  int i,
      j,
      m,
      l,
      n,
      v,
      c;

  double big,
         temp;
  
  theta = dvector(1,NT);
  THETA = dvector(1,NT);
  ERRE = dmatrix(1,NT,1,NR);
  SPH_val = dvector(1,NT*(NR+1)+10);
  init = dvector(1,NT*(NR+1)+10);
  LEG_POL_2c_i = dvector(1,DNMAX*NT);
  fgr = dvector(1,DNMAX*NT*NT*NR*NR);

  for(i=1;i<=NT;i++) {
     theta[i] = (PI/2.0)*((i-1.0)/(NT-1.0));
    }                                            // Definition of angular values for computational grid

  if (argv[1][1]=='f') {sscanf(argv[2],"%s",file_name);}
    else {
      printf("ERROR: please give file with spherical model as input. The computation will be stopped.");
      return 0;
    }
  
  load_sphere(file_name);

  cntr = 0;
  modnum = 0;
  obl=1.0;
      
  for(i=1;i<=NT*(NR+1)+10;i++) init[i] = SPH_val[i];
      

  print_header();
  
  for (i=1;i<=NT*(NR+1)+10;i++) {
      printf("%f ", init[i]);
      if(i==NT*(NR+1)+10) printf("\n");
      }
  
  make_grid(init);

  leg_polyn();

  obl = step;
  modnum=1;

  newt(init,NT*(NR+1)+2,&CKC,vecfunc);

  if(CKC==0) {
    PHYS_QUANT(init,modnum);
    
    for (i=1;i<=NT*(NR+1)+10;i++) {
        printf("%lf ", init[i]);
        if(i==NT*(NR+1)+10) printf("\n");
        }
    }

  else {printf("ERROR IN N-R METHOD\n");
        return 0;}

      
  while(cntr==0 || cntr==6){
            
  obl*=step;
  modnum += 1;

  make_grid(init);

  leg_polyn();
            
  newt(init,NT*(NR+1)+2,&CKC,vecfunc);

  if(CKC==0) {
      PHYS_QUANT(init,modnum);
    
      for (i=1;i<=NT*(NR+1)+10;i++) {
            printf("%lf ", init[i]);
            if(i==NT*(NR+1)+10) printf("\n");
            }
    }

  else {printf("ERROR IN N-R METHOD\n");
        break;}
  
  }

  
  
  return 0;

   
}

