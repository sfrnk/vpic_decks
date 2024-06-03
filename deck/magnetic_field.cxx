// ============================================================================
// Code to compute magnetic field and flux function for thin coils
//
// Some conventions:
// * Coil lies in y-z plane with center axis along x.
// * Positive current flows counter-clockwise about x.
// * Argument lists take coil parameters first, then (x,y,z).
//
// References:
// * J.D.Jackson, Classical Electrodynamics, Chapter 5
// * Simpson+ NASA memo with convenient formulae
//   https://ntrs.nasa.gov/api/citations/20140002333/downloads/20140002333.pdf
//
// --Aaron Tran, started 2023 Dec 08
// ============================================================================

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Special math functions (elliptic integrals) for coil B-fields.
#if __cpp_lib_math_special_functions >= 201603
  // Expected common use case
  #define COMP_ELLINT_1 std::comp_ellint_1
  #define COMP_ELLINT_2 std::comp_ellint_2
#else
  // Older compilers
  #include <tr1/cmath>
  #define COMP_ELLINT_1 std::tr1::comp_ellint_1
  #define COMP_ELLINT_2 std::tr1::comp_ellint_2
#endif


// Returns the magnetic flux function psi(x,y,z) for a thin coil with axial
// position xc, radius rc, and current Ic.
// Code is not tested for vanishing rc->0.
// Units: SI(MKS) with magnetic permeability mu0 implicitly hard-coded to 1.
double psi_coil(
    double xc, double rc, double Ic,
    double x, double y, double z
) {
  double rho = sqrt(y*y + z*z);
  if (rho == 0.) {
    // In this case, asq=bsq --> ksq=0 --> division by zero errors.
    // Use \lim_{\rho->0} A_\phi = 0, so \psi = 0 too.
    return 0.;
  }

  double asq = rc*rc + y*y + z*z + (x-xc)*(x-xc) - 2*rc*rho;
  double bsq = rc*rc + y*y + z*z + (x-xc)*(x-xc) + 2*rc*rho;
  double ksq = 1. - asq/bsq;
  double ellipk = COMP_ELLINT_1(sqrt(ksq));
  double ellipe = COMP_ELLINT_2(sqrt(ksq));

  double Aphi = Ic/M_PI * rc/sqrt(bsq) * (
      (2./ksq - 1.) * ellipk
      - 2./ksq * ellipe
  );
  return rho*Aphi;
}


// Returns the axial magnetic field Bx(x,y,z) for a thin coil with axial
// position xc, radius rc, and current Ic.
// Code is not tested for vanishing rc->0.
// Units: SI(MKS) with magnetic permeability mu0 implicitly hard-coded to 1.
double bx_coil(
    double xc, double rc, double Ic,
    double x, double y, double z
) {
  double rho = sqrt(y*y + z*z);
  double asq = rc*rc + rho*rho + (x-xc)*(x-xc) - 2*rc*rho;
  double bsq = rc*rc + rho*rho + (x-xc)*(x-xc) + 2*rc*rho;
  double ksq = 1. - asq/bsq;
  double ellipk = COMP_ELLINT_1(sqrt(ksq));
  double ellipe = COMP_ELLINT_2(sqrt(ksq));

  // Numerator differs for axial vs radial fields
  double numerator = asq*ellipk + (rc*rc - (x-xc)*(x-xc) - rho*rho)*ellipe;
  double Bx = Ic/M_PI * numerator/(2*asq*sqrt(bsq));
  return Bx;
}


// Returns the radial magnetic field component By(x,y,z) for a thin coil with
// axial position xc, radius rc, and current Ic.
// Code is not tested for vanishing rc->0.
// Units: SI(MKS) with magnetic permeability mu0 implicitly hard-coded to 1.
double by_coil(
    double xc, double rc, double Ic,
    double x, double y, double z
) {
  double rho = sqrt(y*y + z*z);
  if (rho == 0.) {
    // Avoid division by zero; non-axial field is zero by symmetry.
    return 0.;
  }

  double asq = rc*rc + rho*rho + (x-xc)*(x-xc) - 2*rc*rho;
  double bsq = rc*rc + rho*rho + (x-xc)*(x-xc) + 2*rc*rho;
  double ksq = 1. - asq/bsq;
  double ellipk = COMP_ELLINT_1(sqrt(ksq));
  double ellipe = COMP_ELLINT_2(sqrt(ksq));

  // Numerator differs for axial vs radial fields
  double numerator = -asq*ellipk + (rc*rc + (x-xc)*(x-xc) + rho*rho)*ellipe;
  double By = Ic/M_PI * (x-xc)/rho * y/rho * numerator/(2*asq*sqrt(bsq));
  return By;
}


// Returns the radial magnetic field component Bz(x,y,z) for a thin coil with
// axial position xc, radius rc, and current Ic.
// Code is not tested for vanishing rc->0.
// Units: SI(MKS) with magnetic permeability mu0 implicitly hard-coded to 1.
double bz_coil(
    double xc, double rc, double Ic,
    double x, double y, double z
) {
  double rho = sqrt(y*y + z*z);
  if (rho == 0.) {
    // Avoid division by zero; non-axial field is zero by symmetry.
    return 0.;
  }

  double asq = rc*rc + rho*rho + (x-xc)*(x-xc) - 2*rc*rho;
  double bsq = rc*rc + rho*rho + (x-xc)*(x-xc) + 2*rc*rho;
  double ksq = 1. - asq/bsq;
  double ellipk = COMP_ELLINT_1(sqrt(ksq));
  double ellipe = COMP_ELLINT_2(sqrt(ksq));

  // Numerator differs for axial vs radial fields
  double numerator = -asq*ellipk + (rc*rc + (x-xc)*(x-xc) + rho*rho)*ellipe;
  double Bz = Ic/M_PI * (x-xc)/rho * z/rho * numerator/(2*asq*sqrt(bsq));
  return Bz;
}
