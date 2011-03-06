// Evaluate elliptic integrals of the 1st kind
//
// See paper equation 7
// returns \Theta_j
//
// See also Abramowitz, Stegun 17.3.34

#include <iostream>
#include <math.h>

#define GSL_DBL_EPSILON 2.2204460492503131e-16 
#define GSL_SQRT_DBL_EPSILON 1.4901161193847656e-08
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_DBL_MAX        1.7976931348623157e+308

double eval_elliptic(double dr, double dz, double rr_j)
{
  // Precompute values
  double m = 4.0 * rr_j / (dz * dz + dr * dr);    // Argument of elliptic
  double factor = 4.0 / sqrt(dz * dz + dr * dr);
  double val = 0;	// return value
  //std::cout << "m = " << std::scientific << m;
  
  /* [Carlson, Numer. Math. 33 (1979) 1, (4.5)] */
  if(m >= 1.0 || m < 0.0) {
    std::cerr << "Computation of elliptic integral failed (Argument  m = " << m << ">= 1 or < 0)." << std::endl;
  }
  else if(m >= 1.0 - GSL_SQRT_DBL_EPSILON) {
    /* [Abramowitz+Stegun, 17.3.33] */
    const double y = 1.0 - m;
    const double a[] = { 1.38629436112, 0.09666344259, 0.03590092383 };
    const double b[] = { 0.5, 0.12498593597, 0.06880248576 };
    const double ta = a[0] + y*(a[1] + y*a[2]);
    const double tb = -log(y) * (b[0] * y*(b[1] + y*b[2]));
    val = ta + tb;
  }
  else {
    /* This was previously computed as,
     * 
     *         return gsl_sf_ellint_RF_e(0.0, 1.0 - k*k, 1.0, mode, result);
     * 
     *       but this underestimated the total error for small k, since the 
     *       argument y=1-k^2 is not exact (there is an absolute error of
     *       GSL_DBL_EPSILON near y=0 due to cancellation in the subtraction).
     *       Taking the singular behavior of -log(y) above gives an error
     *       of 0.5*epsilon/y near y=0. (BJG) */
    
    double y = 1.0 - m;
    // int status = gsl_sf_ellint_RF_e(0.0, y, 1.0, mode, result);
    double x = 0.0;
    double z = 1.0;
    const double lolim = 5.0 * GSL_DBL_MIN;
    const double uplim = 0.2 * GSL_DBL_MAX;
    const double errtol = 0.001;
    
    if(x+y < lolim || x+z < lolim || y+z < lolim) {
      std::cerr << "Computation of elliptic integral failed (illegal value)." << std::endl;
    }
    else if ( x < uplim && y < uplim && z < uplim) { 
      const double c1 = 1.0 / 24.0;
      const double c2 = 3.0 / 44.0;
      const double c3 = 1.0 / 14.0;
      double xn = x;
      double yn = y;
      double zn = z;
      double mu, xndev, yndev, zndev, e2, e3, s;
      while(1) {
	double lamda;
	double xnroot, ynroot, znroot;
	mu = (xn + yn + zn) / 3.0;
	xndev = 2.0 - (mu + xn) / mu;
	yndev = 2.0 - (mu + yn) / mu;
	zndev = 2.0 - (mu + zn) / mu;
	if (fabs(xndev) < errtol && fabs(yndev) < errtol && fabs(zndev) < errtol) break;
	xnroot = sqrt(xn);
	ynroot = sqrt(yn);
	znroot = sqrt(zn);
	lamda = xnroot * (ynroot + znroot) + ynroot * znroot;
	xn = (xn + lamda) * 0.25;
	yn = (yn + lamda) * 0.25;
	zn = (zn + lamda) * 0.25;
      }
      e2 = xndev * yndev - zndev * zndev;
      e3 = xndev * yndev * zndev;
      s = 1.0 + (c1 * e2 - 0.1 - c2 * e3) * e2 + c3 * e3;
      val = s / sqrt(mu);
      
    }
    else {
      std::cerr << "Computation of elliptic integral failed." << std::endl;
    }
  }
  //std::cout << "\tK(m) computed: " << val << "\t\Theta_j: " << val * factor <<std::endl;
  return val*factor;
}

