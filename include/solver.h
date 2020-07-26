#if !defined(SOLVER_H)
#define SOLVER_H

#include <fftw3.h>

typedef struct {
  // number of grid point
  int nr, nt;
  fftw_plan fwrd, bwrd;
  // right-hand-side of poisson equation
  double *rhs;
  // working arrays, theta-aligned / radial-aligned
  fftw_complex *rhsct, *rhscr;
  // eigenvalues
  double *eigen;
  // radius positions
  double *rs;
  // coefficients of tri-diagonal matrix (dl, d, du)
  double *dl, *d, *du;
  // auxiliary vector
  double *bt;
} solver_t;

extern int init_solver(solver_t **solver, const int nr, const int nt, const double rm, const double rp, const double tm, const double tp);
extern int solve(solver_t *solver);
extern int finalize_solver(solver_t **solver);

#endif // SOLVER_H
