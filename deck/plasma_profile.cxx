/* 
 * Code to tune input plasma profiles for particle init and injection.
 * --Aaron Tran with a lot of help from Doug Endrizzi, started 2023 Dec 01
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Developer notes:
// * Cannot access user_global because it's private in vpic_simulation class.
// * cannot add function to class vpic_simulation post facto.
// Therefore, make your functions pure (avoid dependencies and side effects),
// which is good programming practice anyways.
//
// Remember: code is configuration.

// Warning: the below code does not pay close attention to open vs. closed
// intervals for particle sampling.  Just be careful to not put particles on
// simulation boundaries when choosing your particle injection or
// initialization volumes.

// ============================================================================
// Sampling functions
// ============================================================================

// Uniform random number on (low,high) (open interval)
// copy-pasted from src/vpic/vpic.h to be accessible to functions
// declared outside of vpic_simulation class.
inline double uniform( rng_t * rng, double low, double high ) {
  double dx = drand( rng );
  return low*(1-dx) + high*dx;
}


// Rejection sample for a user-input PDF value within [0,1].
// The rng is typically provided by a vpic_simulation class attribute that is
// visible to input deck code (see src/vpic/vpic.h).
bool rejection_sample(rng_t* rng, double pdf_value) {
  double coin = uniform(rng, 0, 1.0);
  return coin < pdf_value;
}


// ============================================================================
// Interpolation functions
// ============================================================================

// Linearly interpolate y(xv) given strictly monotonic sample points (x,y(x)).
// Strict bounds checking is performed, so this function could be less
// efficient in performance-sensitive code.
// Monotonicity of x is NOT checked, code results may not be predictable
// when you have some x values equal or non-monotonic...
// Spot-check tested by ATr 2023dec11, see test code in:
//   ~/wham/proj/20231211_particle_init_using_flux_function/interpolation_test
// on personal computer.
double interp(double xv, std::vector<double>& x, std::vector<double>& y) {

  // Argument checks
  if (x.size() != y.size()) {
    MESSAGE(("User deck interp got unequal sample arrays, aborting."));
    mp_abort(1000);
  }

  // What monotonic behavior to expect for the rest of the array?
  bool ascending;
  if (x[0] < x[1]) {
    ascending = true;
  } else if (x[0] > x[1]) {
    ascending = false;
  }

  double yv;
  bool found = false;
  for (size_t ii{0}; ii < x.size()-1; ++ii) {
    double x0 = x[ii];
    double x1 = x[ii+1];
    double y0 = y[ii];
    double y1 = y[ii+1];
    // Monotonicity checks
    if (x0 == x1) {
      MESSAGE(("User deck interp got not-strictly-monotonic sample points, aborting."));
      mp_abort(1001);
    } else if (ascending && x0 > x1) {
      MESSAGE(("User deck interp got initially ascending but not monotonic, aborting."));
      mp_abort(1002);
    } else if (!ascending && x0 < x1) {
      MESSAGE(("User deck interp got initially descending but not monotonic, aborting."));
      mp_abort(1003);
    }
    // Check w/equality on both sides of interval to get the very first match;
    // divide by zero is precluded by strict monotonicity checks.
    if ((x0 <= xv && xv <= x1) || (x0 >= xv && xv >= x1)) {
      yv = y0 + (y1-y0) * (xv-x0)/(x1-x0);
      found = true;
    }
  }

  if (!found) {
    MESSAGE(("User deck interp got xv outside range, aborting."));
    mp_abort(1004);
  }

  return yv;
}


// ============================================================================
// Windowing functions
// ============================================================================

// Tukey window on x=[-1,1] with range [0,1], zero elsewhere.
// https://en.wikipedia.org/wiki/Window_function#Tukey_window
// alpha=0 gives cosine (Hann) window, alpha=1 gives box window.
double window_tukey(double x, double alpha) {
  double ax = abs(x);
  if (ax < alpha) {
    return 1.;
  } else if (ax < 1.) {
    return 0.5 + 0.5*cos( (ax-alpha)/(1.-alpha) * M_PI );
  } else {
    return 0.;
  }
}


// cos^2(x^2*pi/2) window on x=[-1,1] with range [0,1], zero elsewhere.
// Suggested by Doug Endrizzi as good approx for radial profiles
// of equilibrium plasmas in a mirror B-field geometry.
double window_cos2x2(double x) {
  if (x*x < 1.) {
    double y = cos(x*x*M_PI/2.);
    return y*y;
  } else {
    return 0.;
  }
}


// Truncated gaussian window on x=[-1,1] with range [0,1], zero elsewhere.
// Suggested by Doug Endrizzi as good approx for axial profiles
// of equilibrium plasmas in a mirror B-field geometry.
double window_trgauss(double x, double sigma) {
  if (x*x < 1.) {
    return exp(-x*x/(2*sigma*sigma));
  } else {
    return 0;
  }
}


// ============================================================================
// Volume selection functions
// ============================================================================

// Sharp-edged cylinder
bool within_cylinder(
    double x, double y, double z,
    double Lxp, double Lzp
) {
  return (abs(x) < Lxp && sqrt(y*y+z*z) < Lzp);
}

// Tapered cylinder PDF with embedded rejection sampling scheme
double tapered_cylinder_pdf(
    double x, double y, double z,
    double Lxp, double Lzp,
    double alphax, double alphaz
) {
  double faxial = window_tukey( x/Lxp, alphax );
  double frho   = window_tukey( sqrt(y*y+z*z)/Lzp, alphaz );
  return frho*faxial;
}
