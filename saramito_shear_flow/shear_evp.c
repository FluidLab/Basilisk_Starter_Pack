#include "navier-stokes/centered.h"
#include "two-phase.h" // I didnt need to use two-phase here to be honest
#include "log-conform.h"
#include "fluidlab_pack.h"

// Comment or uncomment the line below to do Oscillatory shear or start-up shear flow
// #define OSCILLATORY 1

// Saramito parameters
double Re = 0.1; // Reynolds number
double Bi = 3.0; // Bingham number
double Wi = 0.1; // Weissenberg number
double beta = 0.1111; // Viscosity ratio
double regularization = 1e-5; // Regularization parameter

scalar lambdav[], mupv[];
scalar yielded[];

/**
The top and bottom boundary conditions are those of a Couette flow. */
#ifdef OSCILLATORY
  u.t[top] = dirichlet (y*cos(t));
  u.t[bottom] = dirichlet (y*cos(t));
#else
  u.t[top] = dirichlet (y);
  u.t[bottom] = dirichlet (y);
#endif

// Time to end the simulation
double final_time;

int main (int argc, char * argv[]) {

  // Reading some parameters from command line
  if( argc<4 ) {
    printf("Provide parameters please...\n");
    exit(0);
  }
  Wi = atof( argv[1] );
  Bi = atof( argv[2] );
  beta = atof( argv[3] );

  // Creating the domain (rectangle [-10, 10] x [-1, 1])
  // I will mask the y-direction in the init event
  size(20.0);
  origin(-10.0, -10.0);
  periodic(right);
  init_grid(256);

  // Fluid 1 == Fluid 2 (confined)
  rho1 = rho2 = Re;
  mu1 = mu2 = beta;

  /**
  The viscoelastic fields will be set below. */
  mup = mupv;
  lambda = lambdav;

  // Making the poisson solver tolerance a little bit tighter
  TOLERANCE = 1e-4;

  #ifdef OSCILLATORY
    OpenSimulationFolder("oscillatory_Wi%g_Bi%g_beta%g_Re%g", Wi, Bi, beta, Re);
    int number_periods = 5.0;
    final_time = number_periods*2.0*M_PI;
  #else
    OpenSimulationFolder("startup_Wi%g_Bi%g_beta%g_Re%g", Wi, Bi, beta, Re);
    final_time = 100.0; // NOTE: You will a need longer final_time if Wi is really large
  #endif

  // Running the simulation
  run();

  // Closing the folder and open files
  CloseSimulationFolder();
}

// === We need smaller timesteps during the transient period to properly resolve the stresses
// === We can make it bigger as we start approaching the steady state
// #ifndef OSCILLATORY
event control_timestep(i++) {
  DT = 0.001;

  // For startup flow we need a very small timestep at early times to capture the transient curve well
  // (I can probably do this in a smarter way, adjusting DT based on Wi as a timescale)
  #ifndef OSCILLATORY
    if( t<1e-02 )
      DT = .0001;
    else if( t<0.1 )
      DT = 0.001;
    else
      DT = 0.01;
  #endif
}
// #endif

scalar txx_old[], txy_old[], tyy_old[];
int cell_index = -1;
event logfile (i += 1; t<=final_time) {

  /// === Getting the polymeric stress at the center of the domain (should be constant everywhere anyway)
  double txx = interpolate(tau_p.x.x, 0.0, 0.0);
  double txy = interpolate(tau_p.x.y, 0.0, 0.0);
  double tyy = interpolate(tau_p.y.y, 0.0, 0.0);

  /// === Checking how much tau_p has changed since the last saved state
  double change_txx = change(tau_p.x.x, txx_old);
  double change_txy = change(tau_p.x.y, txy_old);
  double change_tyy = change(tau_p.y.y, tyy_old);
  double difference = change_txx + change_txy + change_tyy;

  // Printing some stuff to a log file
  printf("Time step (%d, %lf, %g): %e ... %e %e %e\n", i, t, dt, difference, txx, txy, tyy);
  PrintLog("%d %lf %e %e %e %e\n", i, t, difference, txx, txy, tyy);
}

event print_solution(t+=10.0) {
  /// === Printing the properties below into a VTK file. Can be viewed in Paraview
  scalar *list = {u.x, p, tau_p.x.x, tau_p.x.y, tau_p.y.y};
  const char *list_names[] = {"vel-u", "pressure", "taup_xx", "taup_xy", "taup_yy"};
  PrintMeshVTK_Binary_Float(i, t, list, list_names);
}

event init (t = 0) {
  mask(y>=1.0 ? top : none);
  mask(y<=-1.0 ? bottom : none);

  // Initializing volume fractions and also the velocity profile (du/dy = 1)
  foreach() {
    f[] = ( y<1.0 && y>-1.0 ) ? 1.0 : 0.0;
    u.x[] = y;
  }

}

/// === Updating the EVP properties
event properties (i++) {

  foreach() {
    /// === Norm of the polymeric stress tensor
    double trace_taup = tau_p.x.x[] + tau_p.y.y[];
    double tau_dev_xx = tau_p.x.x[] - 0.5*trace_taup; 
    double tau_dev_yy = tau_p.y.y[] - 0.5*trace_taup; 
    double tau_dev_xy = tau_p.x.y[];
    double norm_tau_dev = sqrt( 1.0*(tau_dev_xx*tau_dev_xx + 2.0*tau_dev_xy*tau_dev_xy + tau_dev_yy*tau_dev_yy) );
    double yield_term = 1.0;
    yielded[] = 1.0;
    if( Bi>0.0 ) {
      yield_term = (norm_tau_dev < (Bi + regularization)) ? regularization : (norm_tau_dev - Bi)/norm_tau_dev;
      yielded[] = (norm_tau_dev < (Bi + regularization)) ? 0.0 : 1.0;
    }
    mupv[] = (1. - beta)*clamp(f[],0,1)/(yield_term);
    lambdav[] = (Wi/yield_term)*clamp(f[],0,1);
  }
}


