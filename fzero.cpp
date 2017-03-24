/*
  A fzero routine copied from fzero.m on MATLAB
 */

#include <cmath>
#include <algorithm>
#include <limits>

typedef double (*FzeroObj)(double, const void *);

double fzero(FzeroObj FunFcn, const void *param, double a, double b,
             double tol, int maxiter, int &info) {
  info = 0;
  double fa = FunFcn(a, param);
  double fb = FunFcn(b, param);

  if (fa*fb > 0) {
    // the sign of the initial fa and fb should be different
    info = 1;
    return std::numeric_limits<double>::quiet_NaN();
  }

  double fc = fb;
  double c=0, d=0, e=0, m=0, toler=0, s=0, q=0, p=0, r=0;
  int iter = 0;
  while (fb != 0 and a != b) {
    // Insure that b is the best result so far, a is the previous
    // value of b, and c is on the opposite side of the zero from b.
    if ((fb > 0) == (fc > 0)) {
      c = a;
      fc = fa;
      d = b - a;
      e = d;
    }

    if (std::abs(fc) < std::abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    // convergence test and possible exit
    m = 0.5 * (c - b);
    toler = 2.0 * tol * std::max(std::abs(b), 1.0);
    if ((std::abs(m) <= toler) || (fb == 0.0))
      break;

    // Choose bisection or interpolation
    if ((std::abs(e) < toler) || (std::abs(fa) <= std::abs(fb))) {
      // Bisection
      d = m;
      e = m;
    } else {
      // Interpolation
      s = fb / fa;
      if (a == c) {
	// Linear interpolation
	p = 2.0 * m * s;
	q = 1.0 - s;
      } else {
	// Inverse quadratic interpolation
	q = fa / fc;
	r = fb / fc;
	p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
	q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0)
	q = -q;
      else
	p = -p;
      // Is interpolated point acceptable
      if ((2.0 * p < 3.0 * m * q - std::abs(toler * q)) && (p < std::abs(0.5 * e * q))) {
	e = d;
	d = p / q;
      } else {
	d = m;
	e = m;
      }
    } // Interpolation

    // Next point
    a = b;
    fa = fb;
    if (std::abs(d) > toler)
      b += d;
    else if (b > c)
      b -= toler;
    else
      b += toler;

    fb = FunFcn(b, param);
    iter = iter + 1;
    if (maxiter < iter)
      break;
  } // Main loop

  return b;
}
