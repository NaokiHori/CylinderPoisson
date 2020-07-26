#include "solver.h"

#define NR 32   // number of grids in r direction
#define NT 64   // number of grids in t direction

#define RM 0.5   // r-minus edge of the coordinate
#define RP 1.0   // r-plus  edge of the coordinate

#define TM 0.0                 // theta-minus
#define TP 6.283185307179586   // theta-plus

static int test_init(double *p, double *q) {
  // p: answer
  // q: right-hand-side of poisson equation
  const double dr=(RP-RM)/NR;
  const double dt=(TP-TM)/NT;
  double C1, C2;
  double r, t;
  int i, j;
  // initialize input (q)
  for(i=0;i<NR;i++){
    r=RM+(i+0.5)*dr;
    for(j=0;j<NT;j++){
      t=1.*j*dt;
      q[i*NT+j]=r*r*cos(2*t);
    }
  }
  // initialize answer (q)
  for(i=0;i<NR;i++){
    r=RM+(i+0.5)*dr;
    for(j=0;j<NT;j++){
      t=1.*j*dt;
      C1=1./6.*pow(RM*RP, 3.)/(pow(RP, 4.)-pow(RM, 4.))*(-pow(RP, 3.)/pow(RM, 3.)+pow(RM, 3.)/pow(RP, 3.));
      C2=1./6.*pow(RM*RP, 3.)/(pow(RP, 4.)-pow(RM, 4.))*(-pow(RP, 3.)*RM+RP*pow(RM, 3.));
      p[i*NT+j]=(1./12.*pow(r, 4.)+C1*pow(r, 2.)+C2/pow(r, 2.))*cos(2.*t);
    }
  }
  return 0;
}

static int test_check_answer(const double *answer, const double *result){
  // L2 norm
  int i, j;
  double err=0.;
  for(i=0;i<NR;i++){
    for(j=0;j<NT;j++){
      err+=pow(answer[i*NT+j]-result[i*NT+j], 2.);
    }
  }
  printf("%.1e\n", sqrt(err/NR/NT));
  return 0;
}

int main(void) {
  solver_t *solver = NULL;
  double *p=malloc(sizeof(double)*NR*NT);
  init_solver(&solver, NR, NT, RM, RP, TM, TP);
  test_init(p, solver->rhs);
  solve(solver);
  test_check_answer(p, solver->rhs);
  finalize_solver(&solver);
  free(p);
  return 0;
}

