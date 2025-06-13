#include <time.h>

#include "headers/centered_odd.h"
#define FILTERED
#include "headers/two-phase_odd.h"
#include "tension.h"
#include "headers/vof-tracer-particles.h"
#include "headers/utility_functions.h"

Particles Pin;

/// === Parameters for the odd momentum equation
double Oh_even;
double Oh_odd;
double viscosity_ratio;
double density_ratio;

/// === Controls the mesh refinement levels
int mesh_level = 10;

/// === Variables to measure the energy budget of the simulation
double kinetic_energy = 0.0;
double surface_energy = 0.0;
double dissipated_energy = 0.0;
double total_energy = 0.0;


/// === Wall clock time used per timestep
double cpu_time_used = 0.0;

// Open boundary at bottom
u.n[bottom] = neumann(0);
p[bottom] = dirichlet(0);
pf[bottom] = dirichlet(0);

// Open boundary on right
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

// Open boundary on the left
u.n[left] = neumann(0);
p[left] = dirichlet(0);
pf[left] = dirichlet(0);

// Open boundary at the top
u.n[top] = neumann(0);
p[top] = dirichlet(0);
pf[top] = dirichlet(0);


int main (int argc, char * argv[]) {

  // === Removing the default dimensions from DT and L0, since we are using nondimensional time and length
  L0 = 0.0 [0];
  DT = 0.0 [0];
  
  // === Checking if we received enough input arguments
  if( argc<5 ) {
    printf("Not enough input arguments..\n");
    exit(0);
  }

  /// === Reading parameters from command line arguments
  Oh_even = atof(argv[1]);
  Oh_odd = atof(argv[2]);
  density_ratio = viscosity_ratio = 1e-02;
  TOLERANCE = 1e-04;
  mesh_level = atoi(argv[3]);
  DT = atof(argv[4]);

  /// === Setting the geometry and initial grid size
  size(8.0);
  origin(-4.0, -4.0);
  init_grid(256);

  /// === Fluid 1: Odd droplet
  rho1 = 1.0;
  mu1 = Oh_even;
  mu_odd1 = Oh_odd;

  /// === Fluid 2: Newtonian "gas" fluid outside the droplet. Very small density and viscosity
  rho2 = density_ratio;
  mu2 = Oh_even*viscosity_ratio;
  mu_odd2 = 0.0;

  /// === Setting the surface tension coefficient
  f.sigma = 1.0;

  OpenSimulationFolder(false, "odd_oscillations_DT%g_mesh%d_OhE%g_OhO%g", DT, mesh_level, Oh_even, Oh_odd);

  run();
  CloseSimulationFolder();
}

event init (t = 0) {

  // double pancake_width = 3.0;
  // double pancake_height = M_PI/4.0/pancake_width;
  // fraction (f, min(min(0.5*pancake_width - x, x + 0.5*pancake_width), min(pancake_height - y, y + 0)) );
  // Pin = new_vof_tracer_particles(1, 1); 
  // particles pp = pl[Pin];
  // pp[0].x = 0.0;
  // pp[0].y = 0.001;

  double center_x = 0.0, center_y = 0.0;
  double radius_x = 0.6, radius_y = 0.5;
  fraction (f, - sq(x - center_x)/sq(radius_x) - sq(y - center_y)/sq(radius_y) + 1.0);

  int number_circles = 4, particles_per_circle = 30;
  Pin = new_vof_tracer_particles(number_circles*particles_per_circle, 1); //Assign phase f[] = 1

  double min_particle_radius_factor = 0.1, max_particle_radius_factor = 0.99;
  double delta_particle_radius_factor = (max_particle_radius_factor - min_particle_radius_factor)/(number_circles - 1);
  double delta_particle_theta = 2.0*M_PI/particles_per_circle;
  int index_particle = 0;
  particles pp = pl[Pin];
  for( int i=0; i<number_circles; i++ ) {
    double current_particle_radius_factor = min_particle_radius_factor + i*delta_particle_radius_factor;

    for( int j=0; j<particles_per_circle; j++ ) {
      double current_particle_theta = 0.0 + j*delta_particle_theta;
      pp[index_particle].x = current_particle_radius_factor*radius_x*cos(current_particle_theta) + center_x;
      pp[index_particle].y = current_particle_radius_factor*radius_y*sin(current_particle_theta) + center_y;
      index_particle++;
    }
  }

  particle_boundary (Pin);

  

  // Making sure boundary conditions are updated
  boundary(all);
}


event adapt (i++) {

  // Strong discretization inside the droplet
  scalar noise_inside[];
  foreach()
    noise_inside[] = f[]*noise();
  adapt_wavelet ({noise_inside, u.x, u.y}, (double[]){1e-4, 1e-4, 1e-4}, maxlevel = mesh_level, minlevel = 5);

  boundary(all);
}

event update_energy_budget(i++) {
  // calculate_energy_budget(&kinetic_energy, &surface_energy, &dissipated_energy,
  //                         Re_even*We, Re_even, We);
  // total_energy = kinetic_energy + surface_energy + dissipated_energy;
}

// event logfile (t+=0.001; t<=50.0) {
event logfile (i++; t<=50.0) {
  scalar w_vort[];
  vorticity(u, w_vort);

  vector force_odd[];
  foreach() {
    double dvdx2 = (u.y[-1, 0] - 2.0*u.y[0, 0] + u.y[1, 0])/sq(Delta);
    double dvdy2 = (u.y[0, -1] - 2.0*u.y[0, 0] + u.y[0, 1])/sq(Delta);
    double dudx2 = (u.x[-1, 0] - 2.0*u.x[0, 0] + u.x[1, 0])/sq(Delta);
    double dudy2 = (u.x[0, -1] - 2.0*u.x[0, 0] + u.x[0, 1])/sq(Delta);
    double center_mu_odd = clamp(f[], 0.0, 1.0)*(mu_odd1 - mu_odd2) + mu_odd2;
    force_odd.x[] = - center_mu_odd*(dvdx2 + dvdy2);
    force_odd.y[] = center_mu_odd*(dudx2 + dudy2);
  }

  double average_vorticity = 0.0;
  double average_force_odd_x = 0.0, average_force_odd_y = 0.0;
  
  
  // We begin by calculating the center of mass (xcm, ycm)
  double wt = 0., xcm = 0., ycm = 0.;
  foreach(reduction(+:wt) reduction(+:xcm) reduction(+:ycm) 
          reduction(+:average_vorticity) reduction(+:average_force_odd_x) reduction(+:average_force_odd_y)) {
    wt += f[]*sq(Delta);
    xcm += f[]*x*sq(Delta);
    ycm += f[]*y*sq(Delta);
    average_vorticity += f[]*w_vort[]*sq(Delta);
    average_force_odd_x += f[]*force_odd.x[]*sq(Delta);
    average_force_odd_y += f[]*force_odd.y[]*sq(Delta);
  }
  xcm /= wt;
  ycm /= wt;
  average_vorticity /= wt;
  average_force_odd_x /= wt;
  average_force_odd_y /= wt;

  PrintLog(true, "%d %lf %e %lf %.6f %.6f %.6f %e %e %e %e\n", 
            i, t, dt, cpu_time_used, kinetic_energy, surface_energy, dissipated_energy, total_energy, 
            average_vorticity, average_force_odd_x, average_force_odd_y);
}

/// === Prints mesh and interface into VTK files
int vtk_step = 0;
event print_solution_mesh(t+=0.05) {
// event print_solution_mesh(i++) {

  // Calculating vorticity
  scalar w_vort[];
  vorticity(u, w_vort);

  // Looping to calculate the flow type parameter
  scalar flow_type[];
  vector force_odd[];
  foreach() {
    // Calculating velocity derivatives
    double dudx = (u.x[1, 0] - u.x[-1, 0])/(2.0*Delta);
    double dudy = (u.x[0, 1] - u.x[0, -1])/(2.0*Delta);
    double dvdx = (u.y[1, 0] - u.y[-1, 0])/(2.0*Delta);
    double dvdy = (u.y[0, 1] - u.y[0, -1])/(2.0*Delta);

    // Calculating strain rate tensor
    double D11 = dudx;
    double D22 = dvdy;
    double D12 = 0.5*( dudy + dvdx );
    double norm_strain = sqrt( sq(D11) + sq(D22) + 2.0*sq(D12) );

    // Calculating vorticity tensor
    double O11 = 0.0, O22 = 0.0;
    double O12 = 0.5*( dvdx - dudy );
    double O21 = 0.5*( dudy - dvdx );
    double norm_vorticity = sqrt( sq(O11) + sq(O22) + sq(O12) + sq(O21) );

    flow_type[] = (norm_strain - norm_vorticity)/( norm_strain + norm_vorticity + 1e-12 );

    double dvdx2 = (u.y[-1, 0] - 2.0*u.y[0, 0] + u.y[1, 0])/sq(Delta);
    double dvdy2 = (u.y[0, -1] - 2.0*u.y[0, 0] + u.y[0, 1])/sq(Delta);
    double dudx2 = (u.x[-1, 0] - 2.0*u.x[0, 0] + u.x[1, 0])/sq(Delta);
    double dudy2 = (u.x[0, -1] - 2.0*u.x[0, 0] + u.x[0, 1])/sq(Delta);
    // double center_mu_odd = clamp(f[], 0.0, 1.0)*(mu_odd1 - mu_odd2) + mu_odd2;
    double center_mu_odd = clamp(f[], 0.0, 1.0)*(1.0 - mu_odd2) + mu_odd2;
    force_odd.x[] = - center_mu_odd*(dvdx2 + dvdy2);
    force_odd.y[] = center_mu_odd*(dudx2 + dudy2);
  }


  scalar *list_scalars = {p, w_vort, flow_type, f};
  const char *scalar_names[] = {"p", "vorticity", "flow_type", "fractions"};
  vector *list_vectors = {u, force_odd};
  const char *vector_names[] = {"u", "force_odd"};
  PrintMeshVTK_Binary_Float(vtk_step, t, true, 
                            list_scalar_data=list_scalars, list_scalar_names=scalar_names,
                            list_vector_data=list_vectors, list_vector_names=vector_names);

  scalar *list_all_scalars = {u.x, u.y, p, w_vort, flow_type, force_odd.x, force_odd.y, f};
  const char *all_scalar_names[] = {"ux", "uy", "p", "vorticity", "flow_type", "force_odd_x", "force_odd_y", "fractions"};
  PrintUniformMeshDataDump(vtk_step, t, 100, 100, 
                              list_all_scalars, all_scalar_names,
                              inside_only = true);

  PrintInterfaceVTK(vtk_step, t);

  PrintTracerParticlesVTK(vtk_step, t, Pin);

  vtk_step++;
}


/// === Event for keeping track of run time per timestep
clock_t start_time = -1.0, end_time = -1.0;
event step_cputime(i++) {
  if( pid() )
    return 0;

  end_time = clock();
  cpu_time_used  = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
  start_time = end_time;
}