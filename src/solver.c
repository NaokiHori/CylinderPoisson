#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include "solver.h"


int init_solver(solver_t **solver, const int nr, const int nt, const double rm, const double rp, const double tm, const double tp){
  const double dr=(rp-rm)/nr;
  const double dt=(tp-tm)/nt;
  fftw_iodim dims[1], howmany_dims[1], thowmany_dims[2];
  double r;
  int i, j;
  *solver=(solver_t*)calloc(1, sizeof(solver_t));
  (*solver)->nt=nt;
  (*solver)->nr=nr;
  (*solver)->rhs=(double*)fftw_malloc(nt*nr*sizeof(double));
  (*solver)->rhsct=(fftw_complex*)fftw_malloc((nt/2+1)*nr*sizeof(fftw_complex));
  (*solver)->rhscr=(fftw_complex*)fftw_malloc((nt/2+1)*nr*sizeof(fftw_complex));
  (*solver)->eigen=(double*)calloc(nt/2+1, sizeof(double));
  (*solver)->rs=(double*)calloc(nr, sizeof(double));
  dims[0].n=nt;
  dims[0].is=1;
  dims[0].os=1;
  howmany_dims[0].n=nr;
  howmany_dims[0].is=nt;
  howmany_dims[0].os=nt/2+1;
  (*solver)->fwrd=fftw_plan_guru_dft_r2c(1, dims, 1, howmany_dims, (*solver)->rhs, (*solver)->rhsct, FFTW_ESTIMATE);
  howmany_dims[0].n=nr;
  howmany_dims[0].is=nt/2+1;
  howmany_dims[0].os=nt;
  (*solver)->bwrd=fftw_plan_guru_dft_c2r(1, dims, 1, howmany_dims, (*solver)->rhsct, (*solver)->rhs, FFTW_ESTIMATE);
  for(j=0;j<nt/2+1;j++){
    (*solver)->eigen[j]=-pow(2.*sin(M_PI*j/nt)/dt, 2.);
  }
  (*solver)->dl=(double*)calloc(nr, sizeof(double));
  (*solver)->d =(double*)calloc(nr, sizeof(double));
  (*solver)->du=(double*)calloc(nr, sizeof(double));
  (*solver)->bt=(double*)calloc(nr, sizeof(double));
  for(i=0;i<nr;i++){
    r=rm+0.5*(2*i+1)*dr;
    (*solver)->rs[i]=r;
    (*solver)->dl[i]=r*( (r-0.5*dr)/dr/dr                 );
    (*solver)->d [i]=r*(-(r-0.5*dr)/dr/dr-(r+0.5*dr)/dr/dr);
    (*solver)->du[i]=r*(                  (r+0.5*dr)/dr/dr);
  }
  (*solver)->d[   0]+=(*solver)->dl[   0];
  (*solver)->d[nr-1]+=(*solver)->du[nr-1];
  return 0;
}

static int zgtsv(const int n, const double *dl, const double *d, const double *du, fftw_complex *q, const double eigen, double *bt){
  // solve Ax=q, where A is a tri-diagonal matrix, x the answer, q the right-hand side
  // similar function as zgtsv in Lapack, but assume q=0 when A is singular
  // sub-diagonal, diagonal, and super-diagonal are dl, d+eigen, and du, respectively
  // all these have size n, but dl[0] and du[n-1] are dummy, not used
  double w;
  int i;
  w=1./(d[0]+eigen);
  bt[0]=w*du[0];
  q[0]=w*q[0];
  for(i=1;i<n-1;i++){
    w=1./(d[i]+eigen-dl[i]*bt[i-1]);
    bt[i]=w*du[i];
    q[i]=w*(q[i]-dl[i]*q[i-1]);
  }
  i=n-1;
  w=d[i]+eigen-dl[i]*bt[i-1];
  if(fabs(w)<DBL_EPSILON){
    q[i]=0.;
  }else{
    w=1./w;
    q[i]=w*(q[i]-dl[i]*q[i-1]);
  }
  for(i=n-2;i>=0;i--){
    q[i]-=bt[i]*q[i+1];
  }
  return 0;
}

int solve(solver_t *solver){
  const int nt=solver->nt;
  const int nr=solver->nr;
  double *rhs=solver->rhs;
  fftw_complex *rhsct=solver->rhsct;
  fftw_complex *rhscr=solver->rhscr;
  const double *eigen=solver->eigen;
  const double *rs=solver->rs;
  const double *dl=solver->dl;
  const double *d =solver->d;
  const double *du=solver->du;
  double *bt=solver->bt;
  int i, j;
  for(i=0;i<nr;i++){
    for(j=0;j<nt;j++){
      rhs[i*nt+j]*=rs[i]*rs[i];
    }
  }
  fftw_execute(solver->fwrd);
  for(j=0;j<nt/2+1;j++){
    for(i=0;i<nr;i++){
      rhscr[j*nr+i]=rhsct[i*(nt/2+1)+j];
    }
  }
  for(j=0;j<nt/2+1;j++){
    zgtsv(nr, dl, d, du, rhscr+j*nr, eigen[j], bt);
  }
  for(j=0;j<nt/2+1;j++){
    for(i=0;i<nr;i++){
      rhsct[i*(nt/2+1)+j]=rhscr[j*nr+i];
    }
  }
  fftw_execute(solver->bwrd);
  for(i=0;i<nr;i++){
    for(j=0;j<nt;j++){
      rhs[i*nt+j]/=nt;
    }
  }
  return 0;
}

int finalize_solver(solver_t **solver){
  fftw_destroy_plan((*solver)->fwrd);
  fftw_destroy_plan((*solver)->bwrd);
  fftw_free((*solver)->rhs);
  fftw_free((*solver)->rhsct);
  fftw_free((*solver)->rhscr);
  free((*solver)->eigen);
  free((*solver)->rs);
  free((*solver)->dl);
  free((*solver)->d);
  free((*solver)->du);
  free((*solver)->bt);
  free(*solver);
  return 0;
}

