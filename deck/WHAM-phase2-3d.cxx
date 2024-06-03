/////////////////////////////////////////////////////
//
//   WHAM Phase 2 with z=+/-1.5m and r=+/-0.15m
//   Only partly capture the expanders
//
//////////////////////////////////////////////////////
// Developer notes.
// * Write your comments now. You will not have time to come back and clean it
//   up later. (source: Enzo dev guide)
// * If you change something, change the comments. Now. Wrong comments are
//   worse than no comments. (source: Enzo dev guide)
// * https://google.github.io/styleguide/cppguide.html
// * https://enzo.readthedocs.io/en/latest/developer_guide/ProgrammingGuide.html
//
// History.
//   2023 Sep/Oct - original deck from Ari Le (LANL)
//   2023 Nov - major cleanup by Aaron Tran (UW--Madison)
//   2023 Nov - forked GDT-3d to target WHAM-phase2-3d parameters
//////////////////////////////////////////////////////


// ============================================================================
// User configuration includes and macros
// ============================================================================

#include <vector>

//#include "dump_info.cxx"
#include "injection_aaron.cxx"
#include "magnetic_field.cxx"
#include "plasma_profile.cxx"
//#include "tally_bc.cxx"
//#include "collisions.cxx"
//#include "collisions-lemons.cc" // Disabled on Ari's advice --ATr,2023nov10

// Number of species being injected
#define NSPEC_PARTICLE_INJ 2

// Employ turnstiles to partially serialize the high-volume file writes.
// In this case, the restart dumps.  Set NUM_TURNSTILES to be the desired
// number of simultaneous writes.
#define NUM_TURNSTILES 256


// ============================================================================
// Other macros
// ============================================================================

// Determine which domains are along the boundaries
// using macro from grid/partition.c
# define RANK_TO_INDEX(rank,ix,iy,iz) BEGIN_PRIMITIVE {           \
    int _ix, _iy, _iz;                                            \
    _ix  = (rank);                /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(topology_x);   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(topology_x);   /* ix = ix */                   \
    _iz  = _iy/int(topology_y);   /* iz = iz */                   \
    _iy -= _iz*int(topology_y);   /* iy = iy */                   \
    (ix) = _ix;                                                   \
    (iy) = _iy;                                                   \
    (iz) = _iz;                                                   \
  } END_PRIMITIVE

// Apply field-injection boundary conditions
// using macro from local.c
#define XYZ_LOOP(xl,xh,yl,yh,zl,zh)             \
  for( z=zl; z<=zh; z++ )                       \
    for( y=yl; y<=yh; y++ )                     \
      for( x=xl; x<=xh; x++ )


// ============================================================================
// begin_globals
// ============================================================================
// If you want to use global variables (for example, to store the dump
// intervals for your diagnostics section), it must be done in the globals
// section. Variables declared the globals section will be preserved across
// restart dumps. For example, if the globals section is:
//   begin_globals {
//     double variable;
//   } end_globals
// the double "variable" will be visible to other input deck sections as
// "global->variable". Note: Variables declared in the globals section are set
// to zero before the user's initialization block is executed. Up to 16K
// of global variables can be defined.
// ============================================================================
begin_globals {

  // Intervals for various outputs and checks
  int restart_interval;
  int fields_interval;
  int particle_interval;
  int energies_interval;

  // Restart dumps and simulation termination from input deck
  int rtoggle;
  int quota_check_interval;
  double quota_sec;

  // Collision parameters for Lemons model
  int tstep_coll;
  double cvar;
  double Z_I1;
  double Z_I2;
  double mi_I1;
  double mi_I2;
  double me_coll;
  double nppc_max;

  // Injector and/or open BC model parameters
  bool use_injector;
  int nrestart;
  int nsp;
  double nfac;
  double ubpar;
  double ubperp;
  double Lxp;
  double Lzp;
  double axp;
  double azp;
  double q[NSPEC_PARTICLE_INJ];
  double npcenter[NSPEC_PARTICLE_INJ];
  double relax[NSPEC_PARTICLE_INJ];
  // Persistent storage of injector plasma moments and cell count
  double *ncenter, *ucenter, *pcenter;
  int *ccenter;
  // Boolean flags to track boundary and injector domains
  int left;
  int right;
  int top;
  int bottom;
  int center;
  // User-added globals for particle/field injection
  double vthi;
  double Lx;
  double Ly;
  double Lz;
  double Lzpinit;
  int sort_interval;

  // Variables for new output format
  DumpParameters fdParams;
  DumpParameters hHdParams;
  DumpParameters hBdParams;
  //DumpParameters hedParams;
  std::vector<DumpParameters *> outputParams;

};


// ============================================================================
// begin_initialization
// ============================================================================
// At this point, there is an empty grid and the random number generator is
// seeded with the rank. The grid, materials, species need to be defined.
// Then the initial non-zero fields need to be loaded at time level 0 and the
// particles (position and momentum both) need to be loaded at time level 0.
// ============================================================================
begin_initialization {

  // -------------------------------------------------------------------------
  // WHAM central cell density is 3e13 cm^-3.
  // Deuterium skin depth is 5.88 cm.
  // Mirror throats are +/-98 cm, so 33.33333*di separation.
  //
  // Take x=+/-1.5m and z=+/-0.15m, so Lx = 300cm and Ly = Lz = 30cm.  This
  // omits full expanders, but captures all field lines that thread the magnet
  // warm bores.
  //
  // Injection region:
  //     Plasma radius 10 cm = 1.70*di
  //     NBI width +/-15cm = 2.55*di (see Figure 12 "z_nbi" in Endrizzi paper)
  // -------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // Plasma parameters
  // --------------------------------------------------------------------------

  // Define the system of units for this problem
  double ec   = 1;                // Charge normalization
  double mi   = 1;                // Mass normalization
  double c    = 1;                // Velocity normalization (speed of light)
  double di   = 1;                // Length normalization (ion inertial length)
  double eps0 = 1;                // Permittivity of space

  // This input deck uses Alfvenic units.
  // Choose B0 = 0.86T and n0 = 3e13 cm^-3 in the central cell,
  // where B_code = n_code both are set to 1.  Then,
  // ... central cell wpi/wci = c/vA0 = 123.8.
  // ... vthi/vA0 = 2.857738e-1 translates to T = 10 keV for deuterium ions.

  // Knob to change mirror size in units of di for fixed thermal bounce time.
  // Equivalent to changing plasma beta while fixing the ratio of Larmor radius
  // to plasma dimension: a/rLi = a/di * di/rLi = a/di * 1/sqrt(beta)
  double fLx = 1;

  // Particle initialization parameters
  double vthi = fLx*sqrt(2)*2.857738e-1;  // Ion thermal velocity sqrt(kB*Ti/mi) in units of vA0
  double TiTe = 2*4;                      // Ion-electron temperature ratio

  // HybridVPIC parameters
  double hyb_te = vthi*vthi/TiTe;   // Electron temperature kB*Te/(mi*vA0^2) scaled to ION mass
  double hyb_den      = 1.0;        // Density corresponding to hyb_te (for hyb_gamma!=1)
  double hyb_den_floor_ohm = 0.05;  // Density floor for Ohm's law update
  double hyb_den_floor_pe = 0.05;   // Density floor for electron pressure update
  double hyb_eta      = 1.0e-4;     // Resistivity
  double hyb_hypereta = 1.0e-4;     // Hyper-resistivity
  double hyb_gamma    = 1.0;        // Electron fluid adiabatic index
  double hyb_nsub     = 10.0;       // Number of field subcycles
  int hyb_nsm         = 0;          // Smoothing passes per timestep for ion moments in field_advance/standard/hyb_advance_b.cc
  int hyb_nsmeb       = 0;          // Smoothing passes per timestep for E/B fields in vpic/advance.cc
  int hyb_nsmb        = 0;          // B-field smoothing interval in vpic/advance.cc, 0 disables

  // Other numerical parameters
  double damp      = 0.0;           // Level of radiation damping

  // --------------------------------------------------------------------------
  // Simulation domain parameters
  // --------------------------------------------------------------------------
  double Lx = fLx*50.*di;   // size of box in x dimension
  double Ly = fLx*6.667*di;   // size of box in y dimension
  double Lz = fLx*6.667*di;   // size of box in z dimension

  // Use double datatype (not int) for n{x,y,z} and topology_{x,y,z} to work
  // with FileIO write to "info.bin" for the Fortran "translate.f90"
  // postprocessing script.  The translate script parses "info.bin" as an
  // unformatted binary stream and so cannot currently cope with mixed double
  // and int datatypes. --ATr,2023nov17
  double nx = 384;          // Number of cells in x, y, and z
  double ny = 96;
  double nz = 96;
  double topology_x = 128;  // Number of domains in x, y, and z
  double topology_y = 2;
  double topology_z = 2;

  double hx = Lx/nx;   // cell size in x
  double hy = Ly/ny;   // cell size in y
  double hz = Lz/nz;   // cell size in z

  // --------------------------------------------------------------------------
  // Simulation magnetic geometry
  // --------------------------------------------------------------------------

  // Four coils mimicking HTS+CC(W7A) magnet arrangement with mirror ratio 20.
  // Axial symmetry about x=0 is assumed in multiple parts of the input deck.
  // Hard-coded numbers assume di = 5.88 cm.
  double xhts0 = -fLx*16.667*di;  // axial position +/-98 cm
  double xhts1 = +fLx*16.667*di;
  double rhts = fLx*3.602*di;  // radius 21.17724 cm
  double Ihts = fLx*142.102;

  double xcc0 = -fLx* 3.401*di;  // axial position +/-20 cm
  double xcc1 = +fLx* 3.401*di;
  double rcc = fLx*12.755*di;  // radius 75 cm
  double Icc = Ihts/16.;

  // Set Bcode=1 at cell center, overriding the user's input current value.
  // Only the relative ratio of Icc/Ihts matters.
  double bscale = 2*bx_coil(xhts0,rhts,Ihts,0,0,0) + 2*bx_coil(xcc0,rcc,Icc,0,0,0);
  Ihts /= bscale;
  Icc /= bscale;

  // --------------------------------------------------------------------------
  // Simulation timestepping, duration, and output parameters
  // --------------------------------------------------------------------------

  // Determine the time step
  double dt = fLx*fLx*0.0001;

  // Set runtime and dump interval to nice multiples of
  // tbounce/dt = (Lx/vthi)/dt, where Lx is the half-distance between mirror
  // throats (Endrizzi+ 2023, Section 1).
  // ... tbounce/dt = 16.667 / 0.4041451837 / 0.0001 = 412,401 (this deck)
  // ... tbounce/dt = 90 / (0.02*10) / (0.0612372/10) = 73,485 (Ari's deck)
  double tbounce = 0.5*(xhts1-xhts0) / vthi / dt;
  int intv = (int)(0.1*tbounce);

  // Simulation runtime
  double num_step_user = 80*intv;   // Simulation runtime in steps

  // Wall-clock runtime
  double quota   = 11.85;           // Wall-clock quota in hours
  double quota_sec = quota*3600;    // Wall-clock quota in seconds

  // Output dumps, restart dumps, and other checks
  int status_interval_user  = 400;      // Stdout/stderr timer reports
  int quota_check_interval  = 400;      // Wall-clock runtime quota check
  int restart_interval      = intv;     // Simulation restart dumps
  int fields_interval       = intv;     // Fields/hydro dumps
  int particle_interval     = 0;        // Particle dumps; 0 to disable
  int energies_interval     = 400;      // Scalar diagnostics
  int fields_stride         = 1;        // Stride for field/hydro dumps


  // --------------------------------------------------------------------------
  // Particle sampling and weights
  // --------------------------------------------------------------------------
  double nppc = 200;                      // Average number of macro particle per cell per species
  double sort_interval = 20;              // Sort interval for particles, type double to match src/vpic/vpic.h
  double Npart  = nppc*nx*ny*nz;          // total macro electrons in box
  Npart = trunc_granular(Npart,nproc());  // Make divisible by number of processors
  double qe = -ec*Lx*Ly*Lz/Npart; // Charge per macro electron
  double qi =  ec*Lx*Ly*Lz/Npart; // Charge per macro ion
                                  // For particle deposit "qe" and "qi" serve as WEIGHT factors;
                                  // "ec" also appears in species definition,
                                  // so particle deposit has two factors of
                                  // "ec" which I think is not correct.  But
                                  // "qe" and "qi" propagate elsewhere into
                                  // the code, and ec=1 means that it's not a
                                  // problem, so I am loath to edit this.
                                  // Just beware.  -ATr,2023nov08

  // --------------------------------------------------------------------------
  // Collision parameters for Lemons model
  // --------------------------------------------------------------------------
  // In CGS, variance of tan theta = 2 pi e^4 n_e dt_coll loglambda / (m_ab^2 u^3)
  // Parameters may also be used for Takizuka-Abe model
  int tstep_coll = 5*sort_interval;   // Collision interval (must be even multiple of sort interval)
  double dt_coll = dt*(tstep_coll);   // in (1/wpe), keep tstep_coll < 1
  double nuii_wci = 5.e-4;            // Collision frequency relative to cyclotron frequency
  double cvar = dt_coll*3.0*sqrt(2.0*M_PI)/4.0*pow(vthi,3)*nuii_wci;  // Base collision variance (dimensionless)
  double Z_I1     = 1.0;                // Ion charge
  double Z_I2     = 1.0;
  double mi_I1    = 1.0;              // Ion mass
  double mi_I2    = 1.0;
  double me_coll  = 1.0/2.*1836;      // Electron mass, 1/2 factor for deuterium ions
  double nppc_max   = 20*nppc;        // Max PPC per species for array pre-allocation

  // --------------------------------------------------------------------------
  // Particle initialization parameters
  // --------------------------------------------------------------------------

  // Specify the real-space density distribution for particle initialization
  // See plasma_profile.cxx, magnetic_field.cxx for details.
  // nprofile = 0 selects Tukey-window tapered cylinder
  //             specified by (Lxpinit, Lzpinit, axpinit, azpinit)
  // nprofile = 1 selects a flux-function dependent profile suggested by Doug
  //             Endrizzi, informed by typical equilibrium shapes in mirror
  //             geometries from FBIS/CQL3D runs (Fall 2023).
  // Beware: (psi,B) interpolation is slow, with nlut=10000 it adds ~40 seconds
  // to simulation init time. --ATr,2023dec14
  int nprofile = 1;

  // Tukey-window tapered cylinder, used with nprofile=0
  double Lxpinit = fLx* 12.0 *di;//10.0*di  // Axial half-width
  double Lzpinit = fLx*  1.5 *di;//1.0*di   // Radius
  double axpinit = 0.667;   // Axial taper: 0=box, 1=Hann window
  double azpinit = 0.333;   // Radial taper: 0=box, 1=Hann window

  // Radial cos^2(psi^2) and axial exp[-(B/Bc-1)^2], used with nprofile=1
  // psi=1 is defined at the radius of the warm bore in the HTS magnets
  // Bc = central mag field strength B(x=0,y=0,z=0)
  double rbore = fLx*0.380*di;          // Warm bore radius to select flux surfaces for particle init.
                                        // Some options, all assuming 5.88cm deuterium skin depth and WHAM Phase II geometry.
                                        // 0.353*di = 2.08cm gives central-cell radius 1.58*di =  9.3cm, psi=1.246.
                                        // 0.380*di = 2.23cm gives central-cell radius 1.70*di = 10.0cm, psi=1.444.
                                        // 0.468*di = 2.75cm gives central-cell radius 2.10*di = 12.3cm, psi=2.195.
  double rmsigma = 1./(sqrt(2)*M_PI);   // Stdev of axial density profile in units of scaled mirror ratio minus 1
  constexpr int nlut = 1000;            // Number of look-up table entries for
                                        // (rho,psi,B) flux surface mapping,
                                        // sampled on rho in [0,Ly/2] and rho in [0,rbore]

  // Velocity space control parameters
  bool avoid_loss_cone = true;          // Sample particles only outside loss cone
                                        // If nprofile=0, mirror ratio is approximated using on-axis fields.
                                        // If nprofile=1, mirror ratio is interpolated from (psi,B) look-up table.

  // --------------------------------------------------------------------------
  // Injector and/or open BC model parameters
  // --------------------------------------------------------------------------

  // Enable or disable continuous particle injection
  bool use_injector = false;
  // When restarting simulation with injectors, need to manually specify the
  // restart directory number to the user deck.
  int inj_nrestart = 0;

  // Specify the volume region for particle injection
  // Tukey-window tapered cylinder matched to nprofile=0 particle init
  double Lxp = fLx*  3.06*di;//2.55*di; // Injection volume axial half-width
  double Lzp = fLx*  1.5 *di;//1.0 *di; // ... and radius
  double axp = 0.667;                   // Axial taper: 0=box, 1=Hann window
  double azp = 0.333;                   // Radial taper: 0=box, 1=Hann window

  // Continuous injection parameters
  int nsp = NSPEC_PARTICLE_INJ;     // Number of species being injected
  double nfac = qi/(hx*hy*hz);      // Normalization factor to convert between PPC and density
  double ubpar = sqrt(12.5)*vthi;   // Beam parallel, perp velocity
  double ubperp = sqrt(12.5)*vthi;
  double inj_q[nsp];                // Species charge, with PPC factors baked in...
  inj_q[0] = qi;
  inj_q[1] = qi;
  double inj_npcenter[nsp];         // Target density for center injector, but
  inj_npcenter[0] = 1.0;            // code is hacked to inject both beam+thermal to n0=1
  inj_npcenter[1] = 0;
  double inj_relax[nsp];            // Relaxation factor (time-damping and delay for injection)
  inj_relax[0] = 0.002*nfac;        // Keep nfac to match Ari's original deck,
  inj_relax[1] = 0.002*nfac;        // but the Nppc dependence may not be desirable --ATr,2023nov16

  // --------------------------------------------------------------------------
  // Initialize high level simulation parameters and globals
  // --------------------------------------------------------------------------

  // Globals in vpic/vpic.h that are directly initialized by user
  num_step                      = num_step_user;
  status_interval               = status_interval_user;

  // Intervals for various outputs and checks
  global->restart_interval      = restart_interval;
  global->fields_interval       = fields_interval;
  global->particle_interval     = particle_interval;
  global->energies_interval     = energies_interval; 

  // Restart dumps and simulation termination from input deck
  global->rtoggle               = 0;  // not meant to be set by user, I think --ATr,2023nov16
  global->quota_check_interval  = quota_check_interval;
  global->quota_sec             = quota_sec;

  // Collision parameters for Lemons model
  global->tstep_coll  = tstep_coll;
  global->cvar        = cvar;
  global->Z_I1        = Z_I1;
  global->Z_I2        = Z_I2;
  global->mi_I1       = mi_I1;
  global->mi_I2       = mi_I2;
  global->me_coll     = me_coll;
  global->nppc_max    = nppc_max;

  // Injector and/or open BC model
  global->use_injector = use_injector;
  global->nrestart    = inj_nrestart;
  global->nsp         = nsp;
  global->nfac        = nfac;
  global->ubpar       = ubpar;
  global->ubperp      = ubperp;
  global->Lxp         = Lxp;
  global->Lzp         = Lzp;
  global->axp         = axp;
  global->azp         = azp;
  global->q[0]        = inj_q[0];
  global->q[1]        = inj_q[1];
  global->npcenter[0] = inj_npcenter[0];
  global->npcenter[1] = inj_npcenter[1];
  global->relax[0]    = inj_relax[0];
  global->relax[1]    = inj_relax[1];

  // User-added globals for particle/field injection
  global->vthi          = vthi;
  global->Lx            = Lx;
  global->Ly            = Ly;
  global->Lz            = Lz;
  global->Lzpinit       = Lzpinit;
  global->sort_interval = sort_interval;

  // --------------------------------------------------------------------------
  // Initialize grid
  // --------------------------------------------------------------------------

  // Setup basic grid parameters
  grid->dx = hx;
  grid->dy = hy;
  grid->dz = hz;
  grid->dt = dt;
  grid->cvac = c;
  grid->eps0 = eps0;

  grid->te = hyb_te;
  grid->den = hyb_den;
  grid->den_floor_ohm = hyb_den_floor_ohm;
  grid->den_floor_pe  = hyb_den_floor_pe;
  grid->eta = hyb_eta;
  grid->hypereta = hyb_hypereta;

  grid->gamma = hyb_gamma;
  grid->nsub = hyb_nsub;
  grid->nsm  = hyb_nsm;
  grid->nsmeb = hyb_nsmeb;
  grid->nsmb = hyb_nsmb;


  // Partition a periodic box among the processors sliced uniformly along x,y,z
  define_periodic_grid(-0.5*Lx, -0.5*Ly, -0.5*Lz,            // Low corner
                        0.5*Lx,  0.5*Ly,  0.5*Lz,            // High corner
                        nx, ny, nz,                          // Resolution
                        topology_x, topology_y, topology_z); // Topology

  // Assign Boolean flags to track boundary and injector domains
  int ix, iy, iz, left=0, right=0, top=0, bottom=0, center=0;
  RANK_TO_INDEX( int(rank()), ix, iy, iz );
  if ( ix ==0 ) left=1;
  if ( ix==topology_x-1) right=1;
  if ( iz ==0 ) bottom=1;
  if ( iz ==topology_z-1 ) top=1;

  // "center" specifies ranks that inject particles and does not correspond to
  // a domain boundary, in contrast to the flags "left", "right", etc...
  // This code must follow the periodic grid partitioning for grid->x0,y0,z0
  // to be valid.
  // BEWARE: nx/ny/nz in current scope is global cell count, not local
  for ( int k=1; k<=(grid->nz); k++ ) {
    for ( int j=1; j<=(grid->ny); j++ ) {
      for ( int i=1; i<=(grid->nx); i++ ) {
        // Cell-center position, don't count ghosts.
        double xc = grid->x0 + hx*(i-0.5);
        double yc = grid->y0 + hy*(j-0.5);
        double zc = grid->z0 + hz*(k-0.5);
        if (within_cylinder(xc,yc,zc,Lxp,Lzp)) {
          // TODO WARNING - this has fallen rather out of date with respect to
          // the injector code and needs to be revisited if/when we begin using
          // particle injection --ATr,2023dec14
          center = 1;
        }
      }
    }
  }

  // Boolean flags to track boundary and injector domains
  // must follow the periodic grid partitioning for grid->x0,y0,z0
  // to be valid for "center" assignment.
  global->left = left;
  global->right = right;
  global->top = top;
  global->bottom = bottom;
  global->center = center;

  // Override some of the boundary conditions
  // "pec_fields" is a perfect electrical conductor
  if ( ix==0 )            set_domain_field_bc( BOUNDARY(-1,0,0), pec_fields );
  if ( ix==topology_x-1 ) set_domain_field_bc( BOUNDARY( 1,0,0), pec_fields );
  if ( iy==0 )            set_domain_field_bc( BOUNDARY(0,-1,0), pec_fields );
  if ( iy==topology_y-1 ) set_domain_field_bc( BOUNDARY(0, 1,0), pec_fields );
  if ( iz==0 )            set_domain_field_bc( BOUNDARY(0,0,-1), pec_fields );
  if ( iz==topology_z-1 ) set_domain_field_bc( BOUNDARY(0,0, 1), pec_fields );

  // Absorbing particle boundaries
  // Disable if using tally BCs which get set up later during init
  if ( ix==0 )            set_domain_particle_bc( BOUNDARY(-1,0,0), absorb_particles );
  if ( ix==topology_x-1 ) set_domain_particle_bc( BOUNDARY( 1,0,0), absorb_particles );
  if ( iy==0 )            set_domain_particle_bc( BOUNDARY(0,-1,0), absorb_particles );
  if ( iy==topology_y-1 ) set_domain_particle_bc( BOUNDARY(0, 1,0), absorb_particles );
  if ( iz==0 )            set_domain_particle_bc( BOUNDARY(0,0,-1), absorb_particles );
  if ( iz==topology_z-1 ) set_domain_particle_bc( BOUNDARY(0,0, 1), absorb_particles );

  // --------------------------------------------------------------------------
  // Initialize materials and field arrays
  // --------------------------------------------------------------------------

  // Note: define_material defaults to isotropic materials with mu=1,sigma=0
  // Tensor electronic, magnetic and conductive materials are supported
  // though. See "sample/shapes" for how to define and assign them to regions.
  material_t * vacuum = define_material( "vacuum", 1 );

  // define_material(...) takes arguments
  //     eps{x,y,z} = relative electrical permittivity
  //     mu{x,y,z} = relative magnetic permeability
  //     sigma{x,y,z} = electrical conductivity
  //     zeta{x,y,z} = magnetic conductivity
  // For HybridVPIC, these arguments are (at least) partly overriden!!
  //   first number epsx multiplies E in hyb_advance_e and hyb_hypereta
  //   second number epsy multiplies eta_hyper in hyb_hypereta
  //   third number epsz multiplies eta in hyb_adv_e

  // Coil forces E->0, perfect conductor
  material_t * coil = define_material("coil",0.,0.,1.,
                                             1.,1.,1.,
                                             0.,0.,0.);

  // Layer applies strong damping to incoming E/M waves
  material_t * layer = define_material("layer", 1.,30.,10.,
                                                1.,1.,1.,
                                                0.,0.,0.);

  // If you pass NULL to define field array, the standard field array will
  // be used (if damp is not provided, no radiation damping will be used).
  //
  // Prerequisites: grid and all materials must be defined before the field
  // array can be defined.
  define_field_array( NULL, damp );

  // Materials and boundaries within the domain are implemented using macros
  // from 'deck/wrapper.cc', see also 'sample/shapes'.
  //
  //     Usage: set_region_material( rgn, interior_material, surface_material );
  //     Usage: set_region_bc( rgn, interior, interior_surface, exterior_surface );
  //     Usage: set_region_fields( rgn,
  //                           eqn_ex, eqn_ey, eqn_ez,
  //                           eqn_bx, eqn_by, eqn_bz );
  //
  // The shape "rgn" is an algebraic expression for global position
  // coordinates (x,y,z) which are not declared in the deck, but rather
  // declared and initialized within the macro, which loops over all cells
  // on this rank.
  //
  // All points inside a cell are in a normal region if the cell center is
  // inside the region.

  // Wrap conductor and resistive material around coils
  // Equations specify ellipses centered on HTS coil position
  double xP1 = xhts0;
  double xP2 = xhts1;
  double RP = rhts;
  double LPin  = fLx*2.5*di;
  double LPout = fLx*2.8*di;  //0.25*Lz;
  double eP = 0.5;  // adjust eccentricity of nested ellipses

#define RHO() ( sqrt(y*y + z*z) )
#define R21   ( eP*(x-xP1)*(x-xP1) + (RHO()-RP)*(RHO()-RP) )
#define R22   ( eP*(x-xP2)*(x-xP2) + (RHO()-RP)*(RHO()-RP) )
#define INSIDE_COIL1  (  R21 <  LPin*LPin )
#define INSIDE_COIL2  (  R22 <  LPin*LPin )
#define INSIDE_LAYER1 ( (R21 >= LPin*LPin) && (R21 < LPout*LPout) )
#define INSIDE_LAYER2 ( (R22 >= LPin*LPin) && (R22 < LPout*LPout) )

  set_region_bc(INSIDE_COIL1, reflect_particles, reflect_particles,reflect_particles);
  set_region_bc(INSIDE_COIL2, reflect_particles, reflect_particles,reflect_particles);

  // Experimental code for faster E/B-field advance requires new macro routine
  // for setting region materials.  Internal field_t structs now track
  // epsx,epsy,epsz to avoid indirect memory access during field advance, see:
  //     src/field_advance/field_advance.h
  //     deck/wrapper.cc
  // You, the user, must explicitly initialize all space with the first
  // material defined to set epsx,epsy,epsz everywhere on the grid.
  // --ATr,2023nov19
  hyb_set_region_material(everywhere, vacuum);
  hyb_set_region_material(INSIDE_LAYER1, layer);
  hyb_set_region_material(INSIDE_LAYER2, layer);
  hyb_set_region_material(INSIDE_COIL1, coil);
  hyb_set_region_material(INSIDE_COIL2, coil);
  // Beware: "coil" at domain edges must be declared AFTER "layer" region.
  // I guess we cannot have "layer" at domain edges.
  // The simulation breaks with NaN values in "energies" scalar diagnostic when
  // "layer" is declared AFTER the rho>0.45*Lz material, tested with 720x72x72
  // cells on Purdue Anvil and my personal Mac. --ATr,2023nov18
  hyb_set_region_material(RHO()> 0.45*Lz, coil);
  hyb_set_region_material(RHO()<-0.45*Lz, coil);
  hyb_set_region_material(x<-0.485*Lx, coil);
  hyb_set_region_material(x>0.485*Lx, coil);

#undef RHO
#undef R1
#undef R22
#undef INSIDE_COIL1
#undef INSIDE_COIL2
#undef INSIDE_LAYER1
#undef INSIDE_LAYER2

  // --------------------------------------------------------------------------
  // Load electromagnetic fields
  // --------------------------------------------------------------------------
  // Note: everywhere is a region that encompasses the entire simulation
  // In general, regions are specifed as logical equations (i.e. x>0 && x+y<2)
  //
  // The field macros
  //     set_region_field(...)
  //     set_region_bext(...)
  // take logical expressions phrased using global coordinates (x,y,z)
  // to initialize electric and magnetic fields.
  // --------------------------------------------------------------------------

#define BXHTS(x,y,z) ( bx_coil(xhts0,rhts,Ihts,(x),(y),(z)) + bx_coil(xhts1,rhts,Ihts,(x),(y),(z)) )
#define BYHTS(x,y,z) ( by_coil(xhts0,rhts,Ihts,(x),(y),(z)) + by_coil(xhts1,rhts,Ihts,(x),(y),(z)) )
#define BZHTS(x,y,z) ( bz_coil(xhts0,rhts,Ihts,(x),(y),(z)) + bz_coil(xhts1,rhts,Ihts,(x),(y),(z)) )

#define BXCC(x,y,z)  ( bx_coil(xcc0,rcc,Icc,   (x),(y),(z)) + bx_coil(xcc1,rcc,Icc,   (x),(y),(z)) )
#define BYCC(x,y,z)  ( by_coil(xcc0,rcc,Icc,   (x),(y),(z)) + by_coil(xcc1,rcc,Icc,   (x),(y),(z)) )
#define BZCC(x,y,z)  ( bz_coil(xcc0,rcc,Icc,   (x),(y),(z)) + bz_coil(xcc1,rcc,Icc,   (x),(y),(z)) )

#define PSIHTS(x,y,z) ( psi_coil(xhts0,rhts,Ihts,(x),(y),(z)) + psi_coil(xhts1,rhts,Ihts,(x),(y),(z)) )
#define PSICC(x,y,z)  ( psi_coil(xcc0,rcc,Icc,   (x),(y),(z)) + psi_coil(xcc1,rcc,Icc,   (x),(y),(z)) )

// These macros will be re-used for particle init and info.hdf5 below
#define BX_LOCAL(x,y,z)  ( BXHTS((x),(y),(z)) + BXCC((x),(y),(z))   )
#define BY_LOCAL(x,y,z)  ( BYHTS((x),(y),(z)) + BYCC((x),(y),(z))   )
#define BZ_LOCAL(x,y,z)  ( BZHTS((x),(y),(z)) + BZCC((x),(y),(z))   )
#define PSI_LOCAL(x,y,z) ( PSIHTS((x),(y),(z)) + PSICC((x),(y),(z)) )
#define MAGNITUDE(i,j,k) ( sqrt((i)*(i) + (j)*(j) + (k)*(k)) )
#define B_LOCAL(x,y,z)   ( MAGNITUDE( BX_LOCAL((x),(y),(z)), BY_LOCAL((x),(y),(z)), BZ_LOCAL((x),(y),(z)) ) )

  set_region_field( everywhere, 0, 0, 0,       // Electric field
                                0, 0 ,0 );     // Magnetic field

  // External Magnetic field
  set_region_bext( everywhere, BX_LOCAL(x,y,z), BY_LOCAL(x,y,z), BZ_LOCAL(x,y,z) );

  // --------------------------------------------------------------------------
  // Initialize species
  // --------------------------------------------------------------------------

  double nmax = 2.*Npart/nproc();
  double nmovers = 0.1*nmax;
  double sort_method = 1; //0=in place, 1=out of place

  species_t* ion  = define_species("ion" ,ec,mi,nmax,nmovers,sort_interval,sort_method);
  species_t* beam = define_species("beam",ec,mi,nmax,nmovers,sort_interval,sort_method);
  //species_t* electron = define_species("electron",-ec,me,nmax,nmovers,sort_interval,sort_method);

  // --------------------------------------------------------------------------
  // Load particles
  // --------------------------------------------------------------------------

  // Make random number generators have a different sequence on each node
  //seed_rand( rng_seed*nproc() + rank() );  //Generators desynchronized
  double xmin = grid->x0 , xmax = grid->x0+(grid->dx)*(grid->nx);
  double ymin = grid->y0 , ymax = grid->y0+(grid->dy)*(grid->ny);
  double zmin = grid->z0 , zmax = grid->z0+(grid->dz)*(grid->nz);

  // Outer boundary flux surface sets radial normalization in psi
  const double psiwb = PSI_LOCAL(xhts1,rbore,0);

  // Construct look-up table for B(psi) at central cell and throat.
  //if (rank() == 0) MESSAGE(("******** Particle init (rho,psi,B) mapping start ********"));
  std::vector<double> lut_psic ( nlut );
  std::vector<double> lut_psit ( nlut );
  std::vector<double> lut_bc   ( nlut );
  std::vector<double> lut_bt   ( nlut );
  // Radial profiles of psi and |B|, with linear sampling in y
  // (psi ~ sqrt(y) within coil is not linearly spaced).
  for (int jj=0; jj<nlut; ++jj) {
    double yc = jj * 0.5*Ly/(nlut-1);  // Linearly sample y in [0,0.5*Ly]
    double yt = jj * rbore/(nlut-1);  // Linearly sample y in [0,rbore]
    lut_psic[jj] = PSI_LOCAL(0,     yc, 0);
    lut_psit[jj] = PSI_LOCAL(xhts1, yt, 0);
    lut_bc[jj] = B_LOCAL(0,     yc, 0);
    lut_bt[jj] = B_LOCAL(xhts1, yt, 0);
    //if (rank() == 0) {
    //  MESSAGE(("Central cell (y,psi,B) = (%g, %g, %g) throat (%g, %g, %g)",
    //        yc, lut_psic[jj], lut_bc[jj],
    //        yt, lut_psit[jj], lut_bt[jj]
    //  ));
    //}
  }
  //if (rank() == 0) MESSAGE(("******** Particle init (rho,psi,B) mapping end ********"));

  // Begin rejection sampling of particles
  repeat ( Npart/nproc() ) {

    // Particle position
    double x = uniform(rng(0),xmin,xmax);
    double y = uniform(rng(0),ymin,ymax);
    double z = uniform(rng(0),zmin,zmax);
    // Fields at particle
    double psi = PSI_LOCAL(x,y,z);
    double bx = BX_LOCAL(x,y,z);
    double by = BY_LOCAL(x,y,z);
    double bz = BZ_LOCAL(x,y,z);
    double b = sqrt(bx*bx + by*by + bz*bz);

    // Central and throat fields on particle's flux surface
    // Use on-axis field values as a cheap approximation / default.
    double bc = B_LOCAL(0,0,0);
    double bt = B_LOCAL(xhts1,0,0);
    // Don't use look-up table outside of plasma's bounding flux surface,
    // else linear interpolation will break.
    // Particles with psi/psiwb > 1 should not spawn for nprofile=1 anyways.
    if (nprofile == 1 && abs(psi/psiwb) < 1) {
      bc = interp(psi, lut_psic, lut_bc);
      bt = interp(psi, lut_psit, lut_bt);
    }

    // Density distribution, normalized to f=1 for n=n0
    double f;
    // Tukey-window tapered cylinder
    if (nprofile == 0) {
      f = tapered_cylinder_pdf(x,y,z,Lxpinit,Lzpinit,axpinit,azpinit);
    // Flux-function dependent density profile
    } else if (nprofile == 1) {
      // Don't inject in expanders or beyond warm-bore flux boundary.
      if ((x <= xhts0 || x >= xhts1) || abs(psi/psiwb) >= 1) {
        f = 0;
      } else {
        // Local mirror ratio minus 1, scaled to max possible value on a given
        // flux surface; rmloc in [0,1].
        double rmloc = (b/bc - 1)/(bt/bc - 1);
        f = window_cos2x2(psi/psiwb) * window_trgauss(rmloc,rmsigma);
      }
    }

    // Is candidate particle within inject spatial region?
    if (rejection_sample(rng(0), f)) {

      double ux = normal(rng(0),0,vthi);
      double uy = normal(rng(0),0,vthi);
      double uz = normal(rng(0),0,vthi);

      // Is candidate particle within inject velocity region?
      if (avoid_loss_cone) {
        double uprll = (ux*bx + uy*by + uz*bz)/b;
        double uu = ux*ux + uy*uy + uz*uz;
        // Escape (loss cone) criterion: u^2 > bt/b * uperp^2
        while ( uu > (bt/b) * (uu - uprll*uprll) ) {
          // Resample
          ux = normal(rng(0),0,vthi);
          uy = normal(rng(0),0,vthi);
          uz = normal(rng(0),0,vthi);
          uprll = (ux*bx + uy*by + uz*bz)/b;
          uu = ux*ux + uy*uy + uz*uz;
        }
      }

      inject_particle( ion, x, y, z, ux, uy, uz, qi, 0, 0);
      //inject_particle( electron, x, y, z,
      //                 normal(rng(0),0,vthe),
      //                 normal(rng(0),0,vthe),
      //                 normal(rng(0),0,vthe),-qe, 0, 0);
    }
  }

  // --------------------------------------------------------------------------
  // Initialize internal (particle pairing) collision operator
  // --------------------------------------------------------------------------

  // Disabled in Ari's 3D mirror deck. -ATr,2023nov02

  // // We don't need this correction here for the LPI deck since the reference
  // // density doesn't have the same reduction as in the interface deck. The nppc calcs
  // // should be correct.
  // double cvar0 = cvar/dt_coll;
  // define_collision_op( takizuka_abe("ii", ion,  ion,  entropy, cvar0, tstep_coll ) );
  // define_collision_op( takizuka_abe("ib", ion,  beam, entropy, cvar0, tstep_coll ) );
  // define_collision_op( takizuka_abe("bb", beam, beam, entropy, cvar0, tstep_coll ) );

  // --------------------------------------------------------------------------
  // Initialize absorbing particle boundary handler
  // --------------------------------------------------------------------------

//  double kemax_mult = 20.;  // max ke in tally is kemax_mult * temperature,
//                            // if t_e = 3 keV then kemax = 300 keV
//  particle_bc_t * absorb_tally_bc = define_particle_bc(
//          absorb_tally( species_list, field_array, entropy )
//  );
//  set_absorb_tally_kemax( absorb_tally_bc, ion,  kemax_mult*1.0*vthi*vthi );
//  set_absorb_tally_kemax( absorb_tally_bc, beam, kemax_mult*1.0*vthi*vthi );
//
//  // Set write interval to hydro interval for the time being
//  set_absorb_tally_write_interval( absorb_tally_bc, hydro_interval );
//
//  // Initialize output directories
//  // energy_flux directory also holds accumulated poynting flux
//  if ( rank()==0 ) dump_mkdir("pinhole");
//  if ( rank()==0 ) dump_mkdir("energy_flux");
//
//  // Paint the simulation volume with materials and boundary conditions
//  // all boundaries are i.v.
//#define IV_region ( x<hx || x>Lx-hx || y<-Ly/2+hy || y>Ly/2-hy || z<-Lz/2+hz || z>Lz/2-hz )
//
//  set_region_bc( IV_region, absorb_tally_bc, absorb_tally_bc, absorb_tally_bc);
//
//#undef IV_region

  // --------------------------------------------------------------------------
  // Log diagnostic information about this simulation
  // --------------------------------------------------------------------------

  sim_log ( "***********************************************" );
  sim_log ( "topology_x = " << topology_x );
  sim_log ( "topology_y = " << topology_y );
  sim_log ( "topology_z = " << topology_z );
  sim_log ( "Ti/Te = " << TiTe ) ;
  sim_log ( "num_step = " << num_step );
  sim_log ( "Lx/di = " << Lx/di );
  sim_log ( "Ly/di = " << Ly/di );
  sim_log ( "Lz/di = " << Lz/di );
  sim_log ( "nx = " << nx );
  sim_log ( "ny = " << ny );
  sim_log ( "nz = " << nz );
  sim_log ( "damp = " << damp );
  sim_log ( "nproc = " << nproc () );
  sim_log ( "nppc = " << nppc );
  sim_log ( "di = " << di );
  sim_log ( "Npart = " << Npart );
  sim_log ( "total # of particles = " << 2*Npart );
  sim_log ( "dt = " << dt );
  //sim_log ( "dt*wpe = " << wpe*dt );
  //sim_log ( "dt*wce = " << wce*dt );
  //sim_log ( "dt*wci = " << wci*dt );
  sim_log ( "dx/di = " << Lx/(di*nx) );
  sim_log ( "dy/di = " << Ly/(di*ny) );
  sim_log ( "dz/di = " << Lz/(di*nz) );
  sim_log ( "vthi = " << vthi );
  sim_log ( "nui/wci = " << nuii_wci );
  sim_log ( "nu*dt_coll = " << nuii_wci*dt_coll );
  sim_log ( "energies_interval = " << energies_interval );
  sim_log ( "fields_interval = " << fields_interval );
  sim_log ( "particle_interval = " << particle_interval );
  sim_log ( "restart_interval = " << restart_interval );
  sim_log ( "sort_interval = " << sort_interval );
  sim_log ( "fields_stride = " << fields_stride );
  sim_log ( "***********************************************" );

  // Dump simulation information to file "info"
  if (rank() == 0 ) {
    FILE *fp_info;
    if ( ! (fp_info=fopen("info.dat", "w")) ) ERROR(("Cannot open file."));
    fprintf( fp_info, "***** Simulation parameters *****\n" );
    fprintf( fp_info, "topology_x =           %g\n", topology_x );
    fprintf( fp_info, "topology_y =           %g\n", topology_y );
    fprintf( fp_info, "topology_z =           %g\n", topology_z );
    fprintf( fp_info, "Ti/Te  =               %g\n", TiTe );
    fprintf( fp_info, "num_step =             %i\n", num_step );
    fprintf( fp_info, "Lx/di =                %g\n", Lx/di );
    fprintf( fp_info, "Ly/di =                %g\n", Ly/di );
    fprintf( fp_info, "Lz/di =                %g\n", Lz/di );
    fprintf( fp_info, "nx =                   %g\n", nx );
    fprintf( fp_info, "ny =                   %g\n", ny );
    fprintf( fp_info, "nz =                   %g\n", nz );
    fprintf( fp_info, "damp =                 %g\n", damp );
    fprintf( fp_info, "nproc =                %g\n", nproc() );
    fprintf( fp_info, "nppc =                 %g\n", nppc );
    fprintf( fp_info, "di =                   %g\n", di );
    fprintf( fp_info, "Npart =                %g\n", Npart );
    fprintf( fp_info, "total # of particles = %g\n", 2*Npart );
    fprintf( fp_info, "dt =                   %g\n", dt );
    //fprintf( fp_info, "dt*wpe =               %g\n", wpe*dt );
    //fprintf( fp_info, "dt*wce =               %g\n", wce*dt );
    //fprintf( fp_info, "dt*wci =               %g\n", wci*dt );
    fprintf( fp_info, "dx/di =                %g\n", Lx/(di*nx) );
    fprintf( fp_info, "dy/di =                %g\n", Ly/(di*ny) );
    fprintf( fp_info, "dz/di =                %g\n", Lz/(di*nz) );
    fprintf( fp_info, "vthi =                 %g\n", vthi );
    fprintf( fp_info, "nu/wci =               %g\n", nuii_wci );
    fprintf( fp_info, "nu*dt_coll =           %g\n", nuii_wci*dt_coll );
    fprintf( fp_info, "energies_interval =    %d\n", energies_interval );
    fprintf( fp_info, "fields_interval =      %d\n", fields_interval );
    fprintf( fp_info, "particle_interval =    %d\n", particle_interval );
    fprintf( fp_info, "restart_interval =     %d\n", restart_interval );
    fprintf( fp_info, "sort_interval =        %g\n", sort_interval ); // double not int
    fprintf( fp_info, "fields_stride =        %d\n", fields_stride );
    fprintf( fp_info, "*********************************\n" );
    fclose(fp_info);
  }

  // Dump simulation information to file "info.bin" for translate script
  // which also gets passed to IDL plotting routines, so be careful
  // about changing this too much...
  if (rank() == 0 ) {

    // write binary info file
    FileIO fp_info;
    if ( ! (fp_info.open("info.bin", io_write)==ok) ) ERROR(("Cannot open file."));

    fp_info.write(&topology_x, 1 );
    fp_info.write(&topology_y, 1 );
    fp_info.write(&topology_z, 1 );

    fp_info.write(&Lx, 1 );
    fp_info.write(&Ly, 1 );
    fp_info.write(&Lz, 1 );

    fp_info.write(&nx, 1 );
    fp_info.write(&ny, 1 );
    fp_info.write(&nz, 1 );

    fp_info.write(&dt, 1 );

    fp_info.write(&mi, 1 ); // Translate script expects this number
    fp_info.write(&vthi, 1 );
    fp_info.write(&status_interval, 1 );
    fp_info.close();
  }

  // Dump simulation information to file "info.hdf5"
  if (rank() == 0 ) {
    std::vector<int> ivs;
    std::vector<double> dvs;
    std::vector<char*> inms, dnms;

    // Int parameters
    // ... domain grid and timestepping
    inms.push_back("/num_step"); ivs.push_back( num_step );
    inms.push_back("/nx"); ivs.push_back( nx ); // different from grid->nx !!
    inms.push_back("/ny"); ivs.push_back( ny );
    inms.push_back("/nz"); ivs.push_back( nz );
    inms.push_back("/topology_x"); ivs.push_back( (int)topology_x );
    inms.push_back("/topology_y"); ivs.push_back( (int)topology_y );
    inms.push_back("/topology_z"); ivs.push_back( (int)topology_z );
    // ... dump intervals and stride
    inms.push_back("/fields_interval"); ivs.push_back( fields_interval );
    inms.push_back("/fields_stride"); ivs.push_back( fields_stride );
    inms.push_back("/energies_interval"); ivs.push_back( energies_interval );
    inms.push_back("/particle_interval"); ivs.push_back( particle_interval );
    inms.push_back("/restart_interval"); ivs.push_back( restart_interval );
    inms.push_back("/sort_interval"); ivs.push_back( sort_interval );

    // Double/float parameters
    // ... domain grid and timestepping
    dnms.push_back("/dt"); dvs.push_back( dt );
    dnms.push_back("/dx"); dvs.push_back( grid->dx );
    dnms.push_back("/dy"); dvs.push_back( grid->dy );
    dnms.push_back("/dz"); dvs.push_back( grid->dz );
    dnms.push_back("/Lx"); dvs.push_back( Lx );
    dnms.push_back("/Ly"); dvs.push_back( Ly );
    dnms.push_back("/Lz"); dvs.push_back( Lz );
    dnms.push_back("/x0"); dvs.push_back( grid->x0 ); // rank=0 gives lower left corner for full domain
    dnms.push_back("/y0"); dvs.push_back( grid->y0 );
    dnms.push_back("/z0"); dvs.push_back( grid->z0 );
    // ... magnetic geometry
    dnms.push_back("/bcenter"); dvs.push_back( B_LOCAL(0,0,0) );
    dnms.push_back("/bthroat"); dvs.push_back( B_LOCAL(xhts1,0,0) );
    dnms.push_back("/xhts0"); dvs.push_back( xhts0 );
    dnms.push_back("/xhts1"); dvs.push_back( xhts1 );
    dnms.push_back("/rhts"); dvs.push_back( rhts );
    dnms.push_back("/Ihts"); dvs.push_back( Ihts );
    dnms.push_back("/xcc0"); dvs.push_back( xcc0 );
    dnms.push_back("/xcc1"); dvs.push_back( xcc1 );
    dnms.push_back("/rcc"); dvs.push_back( rcc );
    dnms.push_back("/Icc"); dvs.push_back( Icc );
    // ... particle and fluid properties
    dnms.push_back("/vthi"); dvs.push_back( vthi );
    dnms.push_back("/TiTe"); dvs.push_back( TiTe );  // redundant w.r.t. hyb_te, but convenient...
    dnms.push_back("/hyb_te"); dvs.push_back( hyb_te );
    dnms.push_back("/hyb_den"); dvs.push_back( hyb_den );
    dnms.push_back("/hyb_eta"); dvs.push_back( hyb_eta );
    dnms.push_back("/hyb_hypereta"); dvs.push_back( hyb_hypereta );
    dnms.push_back("/hyb_gamma"); dvs.push_back( hyb_gamma );
    // ... particle initialization and inject
    dnms.push_back("/nppc"); dvs.push_back( nppc );
    dnms.push_back("/ubpar"); dvs.push_back( ubpar );
    dnms.push_back("/ubperp"); dvs.push_back( ubperp );
    dnms.push_back("/use_injector"); dvs.push_back( use_injector );

    // Almost done, write out file to disk
    //dump_info("info.hdf5", ivs, inms, dvs, dnms);
  }

#undef BXHTS
#undef BYHTS
#undef BZHTS
#undef BXCC
#undef BYCC
#undef BZCC
#undef PSIHTS
#undef PSICC
#undef BX_LOCAL
#undef BY_LOCAL
#undef BZ_LOCAL
#undef PSI_LOCAL
#undef MAGNITUDE
#undef B_LOCAL

  /*--------------------------------------------------------------------------
   * New dump definition
   *------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------
   * Set data output format
   *
   * This option allows the user to specify the data format for an output
   * dump.  Legal settings are 'band' and 'band_interleave'.  Band-interleave
   * format is the native storage format for data in VPIC.  For field data,
   * this looks something like:
   *
   *   ex0 ey0 ez0 div_e_err0 cbx0 ... ex1 ey1 ez1 div_e_err1 cbx1 ...
   *
   * Banded data format stores all data of a particular state variable as a
   * contiguous array, and is easier for ParaView to process efficiently.
   * Banded data looks like:
   *
   *   ex0 ex1 ex2 ... exN ey0 ey1 ey2 ...
   *
   *------------------------------------------------------------------------*/

  global->fdParams.format = band;
  sim_log ( "Fields output format = band" );

  global->hHdParams.format = band;
  sim_log ( "Ion species output format = band" );

  global->hBdParams.format = band;
  sim_log ( "Beam species output format = band" );

  //global->hedParams.format = band;
  //sim_log ( "Electron species output format = band" );

  /*--------------------------------------------------------------------------
   * Set stride
   *
   * This option allows data down-sampling at output.  Data are down-sampled
   * in each dimension by the stride specified for that dimension.  For
   * example, to down-sample the x-dimension of the field data by a factor
   * of 2, i.e., half as many data will be output, select:
   *
   *   global->fdParams.stride_x = 2;
   *
   * The following 2-D example shows down-sampling of a 7x7 grid (nx = 7,
   * ny = 7.  With ghost-cell padding the actual extents of the grid are 9x9.
   * Setting the strides in x and y to equal 2 results in an output grid of
   * nx = 4, ny = 4, with actual extents 6x6.
   *
   * G G G G G G G G G
   * G X X X X X X X G
   * G X X X X X X X G         G G G G G G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G   ==>   G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G G G G G G
   * G G G G G G G G G
   *
   * Note that grid extents in each dimension must be evenly divisible by
   * the stride for that dimension:
   *
   *   nx = 150;
   *   global->fdParams.stride_x = 10; // legal -> 150/10 = 15
   *
   *   global->fdParams.stride_x = 8; // illegal!!! -> 150/8 = 18.75
   *------------------------------------------------------------------------*/

  // fields output path and stride
  sprintf(global->fdParams.baseDir, "fields");
  sprintf(global->fdParams.baseFileName, "fields");
  global->fdParams.stride_x       = fields_stride;
  global->fdParams.stride_y       = fields_stride;
  global->fdParams.stride_z       = fields_stride;
  global->outputParams.push_back(&global->fdParams);

  // ion hydro output path and stride
  sprintf(global->hHdParams.baseDir, "hydro");
  sprintf(global->hHdParams.baseFileName, "Hhydro");
  global->hHdParams.stride_x      = fields_stride;
  global->hHdParams.stride_y      = fields_stride;
  global->hHdParams.stride_z      = fields_stride;
  global->outputParams.push_back(&global->hHdParams);

  // beam hydro output path and stride
  sprintf(global->hBdParams.baseDir, "hydro");
  sprintf(global->hBdParams.baseFileName, "Bhydro");
  global->hBdParams.stride_x      = fields_stride;
  global->hBdParams.stride_y      = fields_stride;
  global->hBdParams.stride_z      = fields_stride;
  global->outputParams.push_back(&global->hBdParams);

  // electron hydro output path and stride
  //global->hedParams.baseDir       = "hydro";
  //global->hedParams.baseFileName  = "ehydro";
  //global->hedParams.stride_x      = fields_stride;
  //global->hedParams.stride_y      = fields_stride;
  //global->hedParams.stride_z      = fields_stride;
  //global->outputParams.push_back(&global->hedParams);

  sim_log ( "Fields x-stride " << global->fdParams.stride_x );
  sim_log ( "Fields y-stride " << global->fdParams.stride_y );
  sim_log ( "Fields z-stride " << global->fdParams.stride_z );

  sim_log ( "Ion species x-stride " << global->hHdParams.stride_x );
  sim_log ( "Ion species y-stride " << global->hHdParams.stride_y );
  sim_log ( "Ion species z-stride " << global->hHdParams.stride_z );

  sim_log ( "Beam species x-stride " << global->hBdParams.stride_x );
  sim_log ( "Beam species y-stride " << global->hBdParams.stride_y );
  sim_log ( "Beam species z-stride " << global->hBdParams.stride_z );

  //sim_log ( "Electron species x-stride " << global->hedParams.stride_x );
  //sim_log ( "Electron species y-stride " << global->hedParams.stride_y );
  //sim_log ( "Electron species z-stride " << global->hedParams.stride_z );

  /*--------------------------------------------------------------------------
   * Set output fields
   *
   * It is now possible to select which state-variables are output on a
   * per-dump basis.  Variables are selected by passing an or-list of
   * state-variables by name.  For example, to only output the x-component
   * of the electric field and the y-component of the magnetic field, the
   * user would call output_variables like:
   *
   *   global->fdParams.output_variables( ex | cby );
   *
   * NOTE: OUTPUT VARIABLES ARE ONLY USED FOR THE BANDED FORMAT.  IF THE
   * FORMAT IS BAND-INTERLEAVE, ALL VARIABLES ARE OUTPUT AND CALLS TO
   * 'output_variables' WILL HAVE NO EFFECT.
   *
   * ALSO: DEFAULT OUTPUT IS NONE!  THIS IS DUE TO THE WAY THAT VPIC
   * HANDLES GLOBAL VARIABLES IN THE INPUT DECK AND IS UNAVOIDABLE.
   *
   * For convenience, the output variable 'all' is defined:
   *
   *   global->fdParams.output_variables( all );
   *------------------------------------------------------------------------*/
  /* CUT AND PASTE AS A STARTING POINT
   * REMEMBER TO ADD APPROPRIATE GLOBAL DUMPPARAMETERS VARIABLE

  output_variables( all );

  // Standard VPIC
  output_variables( electric | div_e_err | magnetic | div_b_err |
                    tca      | rhob      | current  | rhof |
                    emat     | nmat      | fmat     | cmat );

  // HybridVPIC
  output_variables(   electric |   magnetic | magnetic0 |
                    smelectric | smmagnetic |      pexx |
                       current |       rhof |      cmat );

  output_variables( current_density  | charge_density |
                    momentum_density | ke_density     | stress_tensor );
  */

  // These have just the most useful things turned on
  //global->fdParams.output_variables( electric | magnetic | current | div_e_err );
  //global->hedParams.output_variables( current_density | charge_density | ke_density | stress_tensor);
  //global->hHdParams.output_variables( current_density | charge_density | ke_density | stress_tensor);

  // Just EM fields for now
  global->fdParams.output_variables( electric | magnetic | magnetic0 );

  // Here we are dumping everthing
  //global->hedParams.output_variables( all );
  global->hHdParams.output_variables( all );
  global->hBdParams.output_variables( all );

  /*--------------------------------------------------------------------------
   * Convenience functions for simlog output
   *------------------------------------------------------------------------*/

  char varlist[512];
  create_field_list(varlist, global->fdParams);

  sim_log ( "Fields variable list: " << varlist );

  //create_hydro_list(varlist, global->hedParams);

  //sim_log ( "Electron species variable list: " << varlist );

  create_hydro_list(varlist, global->hHdParams);

  sim_log ( "Ion species variable list: " << varlist );

  create_hydro_list(varlist, global->hBdParams);

  sim_log ( "Beam species variable list: " << varlist );

  sim_log("*** Finished with user-specified initialization ***");

  // STANDARD VPIC IS DESCRIBED BELOW, HybridVPIC EVOLUTION LOOP DIFFERS.
  //
  // Upon completion of the initialization, the following occurs:
  // - The synchronization error (tang E, norm B) is computed between domains
  //   and tang E / norm B are synchronized by averaging where discrepancies
  //   are encountered.
  // - The initial divergence error of the magnetic field is computed and
  //   one pass of cleaning is done (for good measure)
  // - The bound charge density necessary to give the simulation an initially
  //   clean divergence e is computed.
  // - The particle momentum is uncentered from u_0 to u_{-1/2}
  // - The user diagnostics are called on the initial state
  // - The physics loop is started
  //
  // The physics loop consists of:
  // - Advance particles from x_0,u_{-1/2} to x_1,u_{1/2}
  // - User particle injection at x_{1-age}, u_{1/2} (use inject_particles)
  // - User current injection (adjust field(x,y,z).jfx, jfy, jfz)
  // - Advance B from B_0 to B_{1/2}
  // - Advance E from E_0 to E_1
  // - User field injection to E_1 (adjust field(x,y,z).ex,ey,ez,cbx,cby,cbz)
  // - Advance B from B_{1/2} to B_1
  // - (periodically) Divergence clean electric field
  // - (periodically) Divergence clean magnetic field
  // - (periodically) Synchronize shared tang e and norm b
  // - Increment the time step
  // - Call user diagnostics
  // - (periodically) Print a status message

} // end of begin_initialization



// ============================================================================
// begin_diagnostics
// ============================================================================
begin_diagnostics {

#if defined(__APPLE__) && defined(__MACH__)
  // Check memory usage on local debugging runs only
  if (step() % status_interval == 0) {
    if (rank() == 0) MESSAGE((" ******** Checking free Memory ********"));
    //https://stackoverflow.com/questions/54287611/measure-the-maximum-memory-usage-during-a-function-call
    //https://stackoverflow.com/questions/10509660/getting-getrusage-to-measure-system-time-in-c
    //https://stackoverflow.com/questions/51291429/linux-getrusage-maxrss-maximum-resident-set-size-not-increasing-with-allocatio
    //https://stackoverflow.com/questions/17779570/printing-long-int-value-in-c
    struct rusage result;
    getrusage(RUSAGE_SELF, &result);
    // Mac OS X maxrss in bytes, Linux in kilobytes; see "man getrusage".
    MESSAGE(("Rank %d maxRSS %ld MB",rank(), result.ru_maxrss/1000/1000));
  }
#endif

  /*--------------------------------------------------------------------------
   * NOTE: YOU CANNOT DIRECTLY USE C FILE DESCRIPTORS OR SYSTEM CALLS ANYMORE
   *
   * To create a new directory, use:
   *
   *   dump_mkdir("full-path-to-directory/directoryname")
   *
   * To open a file, use: FileIO class
   *
   * Example for file creation and use:
   *
   *   // declare file and open for writing
   *   // possible modes are: io_write, io_read, io_append,
   *   // io_read_write, io_write_read, io_append_read
   *   FileIO fileIO;
   *   FileIOStatus status;
   *   status= fileIO.open("full-path-to-file/filename", io_write);
   *
   *   // formatted ASCII  output
   *   fileIO.print("format string", varg1, varg2, ...);
   *
   *   // binary output
   *   // Write n elements from array data to file.
   *   // T is the type, e.g., if T=double
   *   // fileIO.write(double * data, size_t n);
   *   // All basic types are supported.
   *   fileIO.write(T * data, size_t n);
   *
   *   // close file
   *   fileIO.close();
   *------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------
   * Data output directories
   * WARNING: The directory list passed to "global_header" must be
   * consistent with the actual directories where fields and species are
   * output using "field_dump" and "hydro_dump".
   *
   * DIRECTORY PATHS SHOULD BE RELATIVE TO
   * THE LOCATION OF THE GLOBAL HEADER!!!
   * (not sure if this is still required... --ATr,2023nov16)
   *------------------------------------------------------------------------*/

  if(step()==0) {
      dump_mkdir("fields");
      dump_mkdir("hydro");
      //dump_mkdir("rundata");
      dump_mkdir("particles");
      dump_mkdir("restart0");
      dump_mkdir("restart1");  // 1st backup
      dump_mkdir("restart2");  // 2nd backup

      global_header("global", global->outputParams);
      //dump_grid("rundata/grid");
      //dump_materials("rundata/materials");
      //dump_species("rundata/species");
  }

  /*--------------------------------------------------------------------------
   * Normal output dump
   *------------------------------------------------------------------------*/

  const int restart_interval  = global->restart_interval;
  const int fields_interval   = global->fields_interval;
  const int particle_interval = global->particle_interval;
  const int energies_interval = global->energies_interval;

  // Scalar diagnostics
  if (energies_interval > 0 && step() % energies_interval == 0) {
    int append = step() == 0 ? 0 : 1;
    hyb_dump_energies("energies.dat", append);
    hyb_dump_particles_count("particles_count.dat", append);
  }

  if (fields_interval > 0 && step() % fields_interval == 0) {
    field_dump(global->fdParams);
    hydro_dump("ion", global->hHdParams);
    hydro_dump("beam", global->hBdParams);
    //hydro_dump("electron", global->hedParams);
  }

  if (particle_interval > 0 && step() % particle_interval == 0) {

    char subdir[36];

    sprintf(subdir,"particles/T.%d",step());
    dump_mkdir(subdir);

    sprintf(subdir,"particles/T.%d/ion",step());
    dump_particles("ion", subdir);

    sprintf(subdir,"particles/T.%d/beam",step());
    dump_particles("beam", subdir);

    //sprintf(subdir,"particles/T.%d/electron",step());
    //dump_particles("electron", subdir);
  }

  /*--------------------------------------------------------------------------
   * Restart dump
   *------------------------------------------------------------------------*/

  // constants used in DUMP_INJECTORS macro
  const int nsp=global->nsp;
  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;

  if (step() > 0 && step() % restart_interval == 0) {
    if (!global->rtoggle) {
      global->rtoggle = 1;
      BEGIN_TURNSTILE(NUM_TURNSTILES){
        checkpt("restart1/restart", 0);
        if (global->use_injector) {
          DUMP_INJECTORS(1);
        }
      } END_TURNSTILE;
    } else {
      global->rtoggle = 0;
      BEGIN_TURNSTILE(NUM_TURNSTILES){
        checkpt("restart2/restart", 0);
        if (global->use_injector) {
          DUMP_INJECTORS(2);
        }
      } END_TURNSTILE;
    }
  }

  /*--------------------------------------------------------------------------
   * Abort simulation when wall clock time exceeds global->quota_sec
   *------------------------------------------------------------------------*/

  const int quota_check_interval  = global->quota_check_interval;
  const double quota_sec          = global->quota_sec;

  // Synchronize the abort across processors by using mp_elapsed(), which is
  // guaranteed to return the same value for all ranks.
  // Check infrequently to avoid costly ALL_REDUCE in mp_elapsed().

  if (step() > 0 && quota_check_interval > 0 && step() % quota_check_interval == 0) {
    if (uptime() > quota_sec) {
      sim_log("Allowed runtime exceeded for this job.  Terminating....\n");

      BEGIN_TURNSTILE(NUM_TURNSTILES){
        checkpt("restart0/restart",0);
        if (global->use_injector) {
          DUMP_INJECTORS(0);
        }
      } END_TURNSTILE;

      sim_log("Restart dump completed.");
      exit(0); // Exit or abort?
    }
  }

} // end of begin_diagnostics


// ============================================================================
// begin_current_injection
// ============================================================================
begin_current_injection {

  // No current injection for this simulation

}


// ============================================================================
// begin_field_injection
// ============================================================================
begin_field_injection {

  const int nx=grid->nx;
  const int ny=grid->ny;
  const int nz=grid->nz;
  int x,y,z;
  const float x0 = grid->x0, hx = grid->dx, Lx=global->Lx;

  double r = 0.01;

  if (global->left) {
    XYZ_LOOP(0,1,0,ny+1,0,nz+1) field(x,y,z).cbx = (1.-r)*field(x,y,z).cbx;
    XYZ_LOOP(0,1,0,ny+1,0,nz+1) field(x,y,z).cbz = (1.-r)*field(x,y,z).cbz;
  }

  if (global->right) {
    XYZ_LOOP(nx,nx+1,0,ny+1,0,nz+1) field(x,y,z).cbx = (1-r)*field(x,y,z).cbx;
    XYZ_LOOP(nx,nx+1,0,ny+1,0,nz+1) field(x,y,z).cbz = (1-r)*field(x,y,z).cbz;
  }

  // ------------------------------------------------------------------------
  // Radial biasing scheme (start)
  // ------------------------------------------------------------------------
  // Create radial electric field by driving azimuthal magnetic field Btheta
  // at domain edges.
  // Continuation of the "build3_3d" and "build4_3d" experiments performed on
  // Wisconsin's CHTC cluster in late Oct 2023,
  // and "20231121.ari.tite9.posbias" and "20231121.ari.tite9.negbias"
  // experiments performed on Anvil.
  // New changes (2023dec14):
  // * scale ramp-up/down times w.r.t. bounce time,
  // * rescale forcing strength for new code B normalization.
  const float y0 = grid->y0, hy = grid->dy, Ly=global->Ly; 
  const float z0 = grid->z0, hz = grid->dz, Lz=global->Lz; 
  // Drive settings
  bool enable_bias = true;
  float st = step();
  float framp;
  const float t0 = 2*8248;  // Linear ramp-up start // 8248 = 0.02 x (bounce time), see "tbounce" definition above
  const float t1 = 4*8248;  // Linear ramp-up end   // changed from 20231214.WHAM.Lx50.posbias --ATr,2024jan14
  const float t2 =88888888; // Linear ramp-down start
  const float t3 =99999999; // Linear ramp-down end
  const float nu = 0.1;     // Krook-like forcing time
  const float b0 = 0.05;    // Fiducial Btheta strength  // changed from 20231214.WHAM.Lx50.posbias --ATr,2024jan14
  const float ncell = 10;   // Number of cells along axis to apply BC, must be < local nx.  // changed from 20231214.WHAM.Lx50.posbias --ATr,2024jan14
  // Loop variables
  float yglob, zglob, rglob, dglob, bytarget, bztarget;

  // Minimum cell size for y/z slices, both y/z axes scaled to [-1,+1]
  dglob = std::min(hy/(0.5*Ly), hz/(0.5*Lz));

  // Drive settings: apply linear ramp up/down to forcing term
  if (st < t0) {
      framp = 0;
  } else if (st < t1) {
      framp = (st-t0)/(t1-t0);  // linear ramp up from t=t0 to t=t1
  } else if (st < t2) {
      framp = 1;  // constant strength from t=t1 to t=t2
  } else if (st < t3) {
      framp = (t3-st)/(t3-t2); // linear ramp down from t=t2 to t=t3
  } else {
      framp = 0;
  }

  if (enable_bias && global->left) {
    for( z=0; z<=nz+1; z++ ) {
      for( y=0; y<=ny+1; y++ ) {
        for( x=0; x<=ncell; x++ ) { // set axial extent of forcing
          // theta=0 at y-axis and increases CCW about x-axis
          // yglob/rglob = cos(theta)
          // zglob/rglob = sin(theta)
          // Global coordinates scaled to [-1,+1] within bounding coil material
          // at rho > 0.45*Lz (see calls to hyb_set_region_material(...)).
          //
          // CHANGED from 20231214.WHAM.Lx50.posbias to better match plasma
          // radius, assuming that expander mirror ratio is approx same as
          // central-cell mirror ratio.  Previously 0.45*5=2.25di, now lowered
          // to be "Lzpinit" value... --ATr,2024jan14
          //yglob = (y0+(y-0.5)*hy)/(0.45*Ly);
          //zglob = (z0+(z-0.5)*hz)/(0.45*Lz);
          yglob = (y0+(y-0.5)*hy)/(global->Lzpinit);
          zglob = (z0+(z-0.5)*hz)/(global->Lzpinit);
          rglob = sqrt(yglob*yglob + zglob*zglob);
  
          if (rglob >= dglob && rglob <= 1.) {
            // Set twist handedness so left and right walls' axial currents have
            // opposing signs on opposite sides of a field line
            bytarget = -(zglob/rglob) * b0 * ( rglob*rglob - rglob );
            bztarget = +(yglob/rglob) * b0 * ( rglob*rglob - rglob );
          } else {
            bytarget = 0;
            bztarget = 0;
          }
  
          field(x,y,z).cby += -nu*framp*(field(x,y,z).cby - bytarget);
          field(x,y,z).cbz += -nu*framp*(field(x,y,z).cbz - bztarget);
  
        }
      }
    }
  } // endif global->left

  if (enable_bias && global->right) {
    for( z=0; z<=nz+1; z++ ) {
      for( y=0; y<=ny+1; y++ ) {
        for( x=nx-ncell; x<=nx+1; x++ ) { // set axial extent of forcing
          // theta=0 at y-axis and increases CCW about x-axis
          // yglob/r = cos(theta)
          // zglob/r = sin(theta)
          // Global coordinates scaled to [-1,+1] within bounding coil material
          // at rho > 0.45*Lz (see calls to hyb_set_region_material(...)).
          //
          // CHANGED from 20231214.WHAM.Lx50.posbias to better match plasma
          // radius, assuming that expander mirror ratio is approx same as
          // central-cell mirror ratio.  Previously 0.45*5=2.25di, now lowered
          // to be "Lzpinit" value... --ATr,2024jan14
          //yglob = (y0+(y-0.5)*hy)/(0.45*Ly);
          //zglob = (z0+(z-0.5)*hz)/(0.45*Lz);
          yglob = (y0+(y-0.5)*hy)/(global->Lzpinit);
          zglob = (z0+(z-0.5)*hz)/(global->Lzpinit);
          rglob = sqrt(yglob*yglob + zglob*zglob);
  
          if (rglob >= dglob && rglob <= 1.) {
            // Set twist handedness so left and right walls' axial currents have
            // opposing signs on opposite sides of a field line
            bytarget = -1 * -(zglob/rglob) * b0 * ( rglob*rglob - rglob );
            bztarget = -1 * +(yglob/rglob) * b0 * ( rglob*rglob - rglob );
          } else {
            bytarget = 0;
            bztarget = 0;
          }
  
          field(x,y,z).cby += -nu*framp*(field(x,y,z).cby - bytarget);
          field(x,y,z).cbz += -nu*framp*(field(x,y,z).cbz - bztarget);
        }
      }
    }
  } // endif global->right

  // ------------------------------------------------------------------------
  // Radial biasing scheme (end)
  // ------------------------------------------------------------------------

}


// ============================================================================
// begin_particle_collisions
// ============================================================================

begin_particle_collisions {
  // No particle collisions for this simulation
  // Can be provided by other include files
}


// ============================================================================
// begin_particle_injection
// ============================================================================
//
// Central-cell fuelling for mirror problem.
//
// Useful reference implementation, where injectors are placed at domain edges
// to create open boundaries with some inflow:
//
//     sample/reconnection/open-collisional/head/open-collisional
//
// ============================================================================
begin_particle_injection {

  if ( ! global->use_injector) {
    return;
  }

  double x, y, z, rho;
  const int nrestart = global->nrestart;
  const double nfac = global->nfac;
  // injection volume parameters
  const double Lxp = global->Lxp;
  const double Lzp = global->Lzp;
  const double axp = global->axp;
  const double azp = global->azp;
  // nsp,nx,ny,nz must reside in this scope for injection.cxx macros
  const int nsp = global->nsp;
  const int nx = grid->nx;
  const int ny = grid->ny;
  const int nz = grid->nz;
  const double hx = grid->dx;
  const double hy = grid->dy;
  const double hz = grid->dz;

  // ------------------------------------------
  // Initialize the injectors on the first call
  // ------------------------------------------
  // Each injector is a set of particle moment arrays in the y-z plane so
  // that local plasma conditions may inform the injected plasma properties.
  // The moment arrays have shape:
  //     scalar (    nspecies,ny,nz)
  //     vector (  3,nspecies,ny,nz)
  //     tensor (3,3,nspecies,ny,nz)
  // Fortran-style indexing subtly skips over ghost cells.  Arrays are
  // managed using macros from ./injection.cxx and src/util/util_base.h,
  // which gloss over linear-array implementation behind the scenes.
  //
  // Could put
  // ------------------------------------------

  static int initted=0;
  if ( !initted ) {

    if (global->center) DEFINE_INJECTOR(center,ny,nz);
    //if (global->left) DEFINE_INJECTOR(left,  ny,nz);
    //if (global->right)DEFINE_INJECTOR(right, ny,nz);

    // New simulation (as opposed to restart)
    if (step()==0) {

      if (global->center) {
        for ( int n=1; n<=nsp; n++ ) { // species
          for ( int k=1; k<=nz; k++ ) { // grid z-axis
            for ( int j=1; j<=ny; j++ ) { // grid y-axis

              // Density
              ncenter(n,k,j) = npcenter(n)/nfac;
              // Bulk velocity, not currently used
              //ucenter(1,n,k,j) = 0;
              //ucenter(2,n,k,j) = 0;
              //ucenter(3,n,k,j) = 0;
              // Pressure tensor, not currently used
              //pcenter(1,1,n,k,j) = 0.5*ncenter(n,k,j)*global->vthi*global->vthi;
              //pcenter(2,2,n,k,j) = pcenter(1,1,n,k,j);
              //pcenter(3,3,n,k,j) = pcenter(1,1,n,k,j);
              //pcenter(1,2,n,k,j)=0;
              //pcenter(1,3,n,k,j)=0;
              //pcenter(2,1,n,k,j)=0;
              //pcenter(2,3,n,k,j)=0;
              //pcenter(3,1,n,k,j)=0;
              //pcenter(3,2,n,k,j)=0;

              // Number of cells along x within the injection region.
              ccenter(n,k,j) = 0;
              for ( int i=1; i<=nx; i++ ) { // grid x-axis
                // Cell-center position, don't count ghosts.
                double xc = grid->x0 + hx*(i-0.5);
                double yc = grid->y0 + hy*(j-0.5);
                double zc = grid->z0 + hz*(k-0.5);
                if (within_cylinder(xc,yc,zc,Lxp,Lzp)) {
                  // Don't count subcell volumes (i.e. tapered density cells),
                  // because tapering is mostly accounted for via the rejection
                  // sampling scheme.
                  //ccenter(n,k,j) += tapered_cylinder_pdf(xc,yc,zc,Lxp,Lzp,axp,azp);
                  ccenter(n,k,j) += 1;
                }
              } // for i

            } // for j
          } // for k
        } // for species
      } // if center

      //if (global->left) {
      // ...
      //}
      //if (global->right) {
      // ...
      //}

    // Restarting simulation
    } else {

      if (global->center) READ_INJECTOR(center,ny,nz,nrestart);
      //if (global->left) READ_INJECTOR(left,  ny,nz,nrestart);
      //if (global->right)READ_INJECTOR(right, ny,nz,nrestart);

    }

    initted = 1;
  } // End of Initialization

  // ----------------------------------
  // Inject particles at every timestep
  // ----------------------------------
  if (global->center) {

    // Used to construct radially-varying mix of ion/beam in
    // the injection region
    const float Lzp=global->Lzp;

    species_t* ion  = find_species_name( "ion",  species_list );
    species_t* beam = find_species_name( "beam", species_list );

    // No species loop needed here, but initialize dummy index to be consistent
    // with other loops/notation for particle injectors
    int n = 1;

    for ( int k=1; k<=nz; k++ ) {
      for ( int j=1; j<=ny; j++ ) {

        // Inject towards a constant central density npcenter(n) ~ n0*Nppc,
        // multiplied by the number of injection cells along this x line.
        // Assumes the same injection region for both thermal/beam ions, i.e.,
        // ccenter(1,k,j) == ccenter(2,k,j).
        // Beware that npcenter(n) is a macro from injection.cxx!
        //
        // For cells on the injector region edge, count particles in the ENTIRE
        // cell when computing the difference of current and target densities.
        // It's assumed that plasma density in edge cells is roughly uniform
        // for the parts of the cell inside and outside the injector region
        // bounding volume.
        // The number of particles to inject is then computed for the entire
        // edge cell, but only place particles if sampled (x,y,z) lies within
        // the injector volume.
        double ninj = npcenter(n)/nfac - ncenter(n,k,j) - ncenter(n+1,k,j);
        ninj = ninj * ccenter(n,k,j);

        while ( ninj > 0.0 ) {
          ninj--;

          // Candidate particle position in global coordinates
          x = grid->x0 + uniform(rng(0),0,1) * (grid->x1 - grid->x0);
          y = grid->y0 + hy*(j-1) + hy*uniform(rng(0),0,1);
          z = grid->z0 + hz*(k-1) + hz*uniform(rng(0),0,1);
          rho = sqrt(y*y + z*z);

          // Is candidate particle within inject region?
          double f = tapered_cylinder_pdf(x,y,z,Lxp,Lzp,axp,azp);
          if (rejection_sample(rng(0), f)) {
            // Inject 50% beam ions at mirror axis, decreasing with radius
            // to <~10% beam (>~90% thermals) at the edge of injection
            // region, for the case of Lzp=Lz/4.5.
            double coin = uniform(rng(0),0,1.0);
            if (coin > 0.5*exp(-2.0*rho*rho/Lzp/Lzp)) {
              // Add a thermal ion
              double ux = normal(rng(0),0,global->vthi);
              double uy = normal(rng(0),0,global->vthi);
              double uz = normal(rng(0),0,global->vthi);
              inject_particle(ion, x, y, z, ux, uy, uz, q(n), 0, 0);

            } else {
              // Add a beam ion
              double signx = uniform(rng(0),0,1);
              double phi = uniform(rng(0),0,2*M_PI);
              double ux = global->ubpar * (signx < 0.5 ? -1 : +1);
              double uy = global->ubperp * cos(phi);
              double uz = global->ubperp * sin(phi);
              inject_particle(beam, x, y, z, ux, uy, uz, q((n+1)), 0, 0 );

            } // end if coin

          } // end if within region

        } // end repeat

      } // end for grid y-axis loop
    } // end for grid z-axis loop
  } // end center injector

  //if (global->left) {
  // ...
  //}

  //if (global->right) {
  // ...
  //}

  // --------------------------------------------------
  // Update the injector moments at every sort interval
  // --------------------------------------------------

#define icell(i,j,k) INDEX_FORTRAN_3(i,j,k,0,nx+1,0,ny+1,0,nz+1)

  if (step() % global->sort_interval == 0) {

    // Sum all particles along x within the center injector region
    if (global->center) {
      for ( int n=1; n<= nsp; n++ ) {

        species_t* species = find_species_id(n-1,species_list);
        double rn = global->relax[n-1];

        for ( int k=1; k<=nz; k++ ) {
          for ( int j=1; j<=ny; j++ ) {

            int npart = 0;
            int ncell = 0;
            double ncellf = 0;

            for ( int i=1; i<=nx; i++ ) {
              // Cell-center position, don't count ghosts.
              // Local domain's lower-left cell has i=j=k=1.
              double xc = grid->x0 + hx*(i-0.5);
              double yc = grid->y0 + hy*(j-0.5);
              double zc = grid->z0 + hz*(k-0.5);
              if (within_cylinder(xc,yc,zc,Lxp,Lzp)) {
                // Count particles in each voxel (no ghosts), only works
                // immediately after a sort
                // The species->partition array is explained in
                // src/species_advance/species_advance.h
                int nstart = species->partition[icell(i,j,k)];
                int nstop  = species->partition[icell(i,j,k)+1];
                npart += (nstop - nstart);
                ncell += 1;
                ncellf += tapered_cylinder_pdf(xc,yc,zc,Lxp,Lzp,axp,azp);
              }
            } // for i

            // Convert raw counts to mean particles per cell and update
            // according to a relaxation factor.
            if (ncellf > 0) {
              ncenter(n,k,j) = (1.0-rn)*ncenter(n,k,j) + rn*npart/ncellf;
            }
            // Number of injection cells along x, which may evolve in time
            // or vary between different species, in principle.
            // Don't count subcell volumes (i.e. tapered density cells),
            // because tapering is mostly accounted for via the rejection
            // sampling scheme.
            ccenter(n,k,j) = ncell;

          } // for j
        } // for k

      } // for n species
    } // if center

    // For open boundary conditions at domain edges, look at
    //   sample/reconnection/open-collisional/head/open-collisional
    // and Daughton,Scudder,Karimabadi (2006,PoP).

    //if (global->left) {
    // ...
    //}

    //if (global->right) {
    // ...
    //}

  } // if sort interval

} // end of begin_particle_injection
