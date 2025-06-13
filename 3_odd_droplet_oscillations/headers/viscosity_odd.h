#include "poisson.h"

struct Viscosity {
  face vector mu;
  face vector mu_odd;
  scalar rho;
  double dt;
};

#if AXI
# define lambda ((coord){1., 1. + dt/rho[]*(mu.x[] + mu.x[1] + \
					    mu.y[] + mu.y[0,1])/2./sq(y)})
#else // not AXI
# if dimension == 1
#   define lambda ((coord){1.})
# elif dimension == 2
#   define lambda ((coord){1.,1.})
# elif dimension == 3
#   define lambda ((coord){1.,1.,1.})
#endif
#endif

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) face vector mu_odd = p->mu_odd;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif
  
  foreach_level_or_leaf (l) {

    /*  In the original viscosity.h they used a foreach_dimension loop to generalize both velocity directions.
        In the odd-viscosity case we will need to individually do each direction, because the symmetry in the x/y equations is broken.
        Note: I will keep the implementation below specific to 2 dimensions.
    */

    w.x[] = (dt/rho[]*(
            2.*mu.x[1, 0]*u.x[1, 0] + 2.*mu.x[]*u.x[-1, 0]
            + mu.y[0,1]*(u.x[0,1] + (u.y[1,0] + u.y[1,1])/4. - (u.y[-1,0] + u.y[-1,1])/4.)
            - mu.y[]*(- u.x[0,-1] + (u.y[1,-1] + u.y[1,0])/4. - (u.y[-1,-1] + u.y[-1,0])/4.)
            - mu_odd.x[1, 0]*( u.y[1, 0] - u.y[0, 0] + (u.x[1, 1] + u.x[0, 1])/4. - (u.x[1, -1] + u.x[0, -1])/4. )
            + mu_odd.x[]*( u.y[0, 0] - u.y[-1, 0] + (u.x[0, 1] + u.x[-1, 1])/4. - (u.x[0, -1] + u.x[-1, -1])/4. )
            - mu_odd.y[0, 1]*( u.y[0, 1] - u.y[0, 0] + (u.x[-1, 1] + u.x[-1, 0])/4. - (u.x[1, 1] + u.x[1, 0])/4. )
            - mu_odd.y[]*( u.y[0, 0] - u.y[0, -1] + (u.x[1, 0] + u.x[1, -1])/4. - (u.x[-1, 0] + u.x[-1, -1])/4. )
            ) 
            + r.x[]*sq(Delta)
            )/(sq(Delta) + dt/rho[]*(2.*mu.x[1, 0] + 2.*mu.x[] + mu.y[0,1] + mu.y[]));

    w.y[] = (dt/rho[]*(
            2.*mu.y[0, 1]*u.y[0, 1] + 2.*mu.y[]*u.y[0, -1]
            + mu.x[1, 0]*(u.y[1, 0] + (u.x[0, 1] + u.x[1,1])/4. - (u.x[0, -1] + u.x[1, -1])/4.)
            - mu.x[]*(- u.y[-1, 0] + (u.x[-1, 1] + u.x[0, 1])/4. - (u.x[-1,-1] + u.x[0, -1])/4.)
            + mu_odd.x[1, 0]*( u.x[1, 0] - u.x[0, 0] - (u.y[1, 1] + u.y[0, 1])/4. + (u.y[1, -1] + u.y[0, -1])/4. )
            - mu_odd.x[]*( u.x[0, 0] - u.x[-1, 0] - (u.y[0, 1] + u.y[-1, 1])/4. + (u.y[0, -1] + u.y[-1, -1])/4. )
            + mu_odd.y[0, 1]*( u.x[0, 1] - u.x[0, 0] - (u.y[-1, 1] + u.y[-1, 0])/4. + (u.y[1, 1] + u.y[1, 0])/4. )
            + mu_odd.y[]*( -(u.x[0, 0] - u.x[0, -1]) - (u.y[1, 0] + u.y[1, -1])/4. + (u.y[-1, 0] + u.y[-1, -1])/4. )
            ) 
            + r.y[]*sq(Delta)
            )/(sq(Delta) + dt/rho[]*(2.*mu.y[0, 1] + 2.*mu.y[] + mu.x[1,0] + mu.x[]));      
      
  }

#if JACOBI
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = (u.x[] + 2.*w.x[])/3.;
#endif
  
#if TRASH
  vector u1[];
  foreach_level_or_leaf (l)
    foreach_dimension()
      u1.x[] = u.x[];
  trash ({u});
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = u1.x[];
#endif
}

// #define ORIGINAL

#ifndef ORIGINAL
static double residual_viscosity (scalar * a, scalar * b, scalar * resl, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) face vector mu_odd = p->mu_odd;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;
  
  /* conservative coarse/fine discretisation (2nd order) */

  /**
  We manually apply boundary conditions, so that all components are
  treated simultaneously. Otherwise (automatic) BCs would be applied
  component by component before each foreach_face() loop. */
  
  boundary ({u});
  
  /*  In the original viscosity.h they used a foreach_dimension loop to generalize both velocity directions.
      In the odd-viscosity case we will need to individually do each direction, because the symmetry in the x/y equations is broken.
      Note: I will keep the implementation below specific to 2 dimensions.
  */

  /// === Begin by calculating the two stress components for the x-direction
  face vector taux[];
  foreach_face(x) { /// Calculating the xx stress component
    taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1, 0])/Delta; 
    taux.x[] += mu_odd.x[]*( - (u.y[] - u.y[-1, 0]) - (u.x[-1, 1] + u.x[0, 1])/4.0 + (u.x[-1, -1] + u.x[0, -1])/4.0 )/Delta; 
  }
  foreach_face(y) { /// Calculating the xy stress component
    taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] + (u.y[1,-1] + u.y[1,0])/4. - (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
    taux.y[] += mu_odd.y[]*( - (u.y[0, 0] - u.y[0, -1]) + (u.x[1,-1] + u.x[1,0])/4. - (u.x[-1,-1] + u.x[-1,0])/4. )/Delta; 
  }

  

  /// === Now calculating the two stress components for the y-direction
  face vector tauy[];
  foreach_face(x) { /// Calculating the yx stress component
    tauy.x[] = mu.x[]*(u.y[] - u.y[-1, 0] + (u.x[-1, 1] + u.x[0, 1])/4. - (u.x[-1,-1] + u.x[0, -1])/4.)/Delta;
    tauy.x[] += mu_odd.x[]*(u.x[] - u.x[-1, 0] - (u.y[-1, 1] + u.y[0, 1])/4. + (u.y[-1,-1] + u.y[0, -1])/4.)/Delta;
  }
  foreach_face(y) { /// Calculating the yy stress component
    tauy.y[] = 2.*mu.y[]*(u.y[] - u.y[0, -1])/Delta;
    tauy.y[] += mu_odd.y[]*( u.x[] - u.x[0, -1] + (u.y[1, 0] + u.y[1, -1])/4.0 - (u.y[-1, 0] + u.y[-1, -1])/4.0 )/Delta; 
  }

  /// Calculating div(tau) and then the residual
  foreach( reduction(max:maxres) ) {
    double div_x = taux.x[1, 0] - taux.x[] + taux.y[0, 1] - taux.y[];
    double div_y = tauy.x[1, 0] - tauy.x[] + tauy.y[0, 1] - tauy.y[];
    res.x[] = r.x[] - lambda.x*u.x[] + dt/rho[]*div_x/Delta;
    res.y[] = r.y[] - lambda.y*u.y[] + dt/rho[]*div_y/Delta;

    if (fabs (res.x[]) > maxres)
      maxres = fabs (res.x[]);
    if (fabs (res.y[]) > maxres)
      maxres = fabs (res.y[]);
  }

  // printf("odd maxres %e\n", maxres);
  // exit(0);
  // getchar();

  return maxres;
}
#else
static double residual_viscosity (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */

  /**
  We manually apply boundary conditions, so that all components are
  treated simultaneously. Otherwise (automatic) BCs would be applied
  component by component before each foreach_face() loop. */
  
  boundary ({u});
  
  foreach_dimension() {
    face vector taux[];
    foreach_face(x)
      taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
    #if dimension > 1
      foreach_face(y)
	taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] + 
			   (u.y[1,-1] + u.y[1,0])/4. -
			   (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
    #endif
    #if dimension > 2
      foreach_face(z)
	taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] + 
			   (u.z[1,0,-1] + u.z[1,0,0])/4. -
			   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
    #endif
    foreach (reduction(max:maxres)) {
      double d = 0.;
      foreach_dimension()
	d += taux.x[1] - taux.x[];
      res.x[] = r.x[] - lambda.x*u.x[] + dt/rho[]*d/Delta;
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres))
    foreach_dimension() {
      res.x[] = r.x[] - lambda.x*u.x[] +
        dt/rho[]*(2.*mu.x[1,0]*(u.x[1] - u.x[])
		  - 2.*mu.x[]*(u.x[] - u.x[-1])
        #if dimension > 1
		  + mu.y[0,1]*(u.x[0,1] - u.x[] +
			       (u.y[1,0] + u.y[1,1])/4. -
			       (u.y[-1,0] + u.y[-1,1])/4.)
		  - mu.y[]*(u.x[] - u.x[0,-1] +
			    (u.y[1,-1] + u.y[1,0])/4. -
			    (u.y[-1,-1] + u.y[-1,0])/4.)
	#endif
        #if dimension > 2
		  + mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
				 (u.z[1,0,0] + u.z[1,0,1])/4. -
				 (u.z[-1,0,0] + u.z[-1,0,1])/4.)
		  - mu.z[]*(u.x[] - u.x[0,0,-1] +
			    (u.z[1,0,-1] + u.z[1,0,0])/4. -
			    (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
	#endif
		  )/sq(Delta);
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
#endif

  // printf("original maxres: %e\n", maxres);
  return maxres;
}
#endif

#undef lambda

trace
mgstats viscosity (vector u, face vector mu, face vector mu_odd, scalar rho, double dt,
		   int nrelax = 4, scalar * res = NULL)
{
  vector r[];
  foreach()
    foreach_dimension()
      r.x[] = u.x[];

  restriction ({mu, mu_odd, rho});
  struct Viscosity p = { mu, mu_odd, rho, dt };
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, nrelax, res);
}

trace
mgstats viscosity_explicit (vector u, face vector mu, face vector mu_odd, scalar rho, double dt)
{
  vector r[];
  mgstats mg = {0};
  struct Viscosity p = { mu, mu_odd, rho, dt };
  mg.resb = residual_viscosity ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  foreach()
    foreach_dimension()
      u.x[] += r.x[];
  return mg;
}
