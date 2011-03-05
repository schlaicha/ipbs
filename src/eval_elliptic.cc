// Evaluate elliptic integrals of the 1st kind
//
// See paper equations 5 - 7
// See also Abramowitz, Stegun 17.3.34

#include <iostream>
#include <math.h>

#define GSL_DBL_EPSILON 2.2204460492503131e-16 
#define GSL_SQRT_DBL_EPSILON 1.4901161193847656e-08

double eval_elliptic(double dr, double dz, double rr_j)
{
  // Precompute values
  double m = -4.0 * rr_j / (dz * dz + dr * dr);    // Argument of elliptic
  
  /* [Carlson, Numer. Math. 33 (1979) 1, (4.5)] */
  if(m >= 1.0) {
    std::cerr << "Computation of elliptic integral failed (Argument >= 1)." << std::endl;
  }
  else if(m >= 1.0 - GSL_SQRT_DBL_EPSILON) {
    /* [Abramowitz+Stegun, 17.3.33] */
    const double y = 1.0 - m;
    const double a[] = { 1.38629436112, 0.09666344259, 0.03590092383 };
    const double b[] = { 0.5, 0.12498593597, 0.06880248576 };
    const double ta = a[0] + y*(a[1] + y*a[2]);
    const double tb = -log(y) * (b[0] * y*(b[1] + y*b[2]));
    const double val = ta + tb;
  }
  else {
    /* This was previously computed as,

         return gsl_sf_ellint_RF_e(0.0, 1.0 - k*k, 1.0, mode, result);

       but this underestimated the total error for small k, since the 
       argument y=1-k^2 is not exact (there is an absolute error of
       GSL_DBL_EPSILON near y=0 due to cancellation in the subtraction).
       Taking the singular behavior of -log(y) above gives an error
       of 0.5*epsilon/y near y=0. (BJG) */

    /*double y = 1.0 - k*k;
    int status = gsl_sf_ellint_RF_e(0.0, y, 1.0, mode, result);
    result->err += 0.5 * GSL_DBL_EPSILON / y;
    return status ;*/
  }
  return 0;
}

