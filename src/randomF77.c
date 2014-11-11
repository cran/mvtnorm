/* $Id: randomF77.c 291 2014-11-11 15:55:04Z thothorn $
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include <R.h>
#include <Rmath.h>

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }
/* double F77_SUB(sqrtqchisq)(double *n, double *p) { return(sqrt(qchisq(p[0], n[0], 0, 0))); } */
double F77_SUB(sqrtqchisqint)(int *n, double *p) { return(sqrt(qchisq(p[0], (double) n[0], 0, 0))); }
