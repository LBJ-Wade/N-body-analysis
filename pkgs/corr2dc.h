#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Boundary conditions */
#define CWRAP(x, y) (x >= y ? x - y : x < 0 ? x + y : x) /* Cell/box wrap to [0,boxsize) */
#define RELWRAP(x, y) (x >= y/2 ? x - y : x < -y/2 ? x + y : x) /* Relative wrap to [-boxsize/2, boxsize/2) */

void init_mesh(long *ll, long *hoc, double *p, long np, int nattr, int nho, double blen);

void corr2d(double *xc_2d, double *p1, long np1, double *p2, long np2, int nbins,
        int nhocells, double blen, double dis_f, int weighted, int ispolar);

void estimator1d(double *xc_1d, double *npair, double *p1, long np1, double *p2, long np2,
        int nbins, int nhocells, double blen, double dis_i, double dis_f, int weighted);

void estimator1dhaloterm(double *xc_1d, double *npair, double *p1, long np1, double *p2, long np2,
        int nbins, int nhocells, double blen, double dis_i, double dis_f, int weighted, int haloterm);

void corr1d(double *xc_1d, double *p1, long np1, double *p2, long np2, int nbins,
        int nhocells, double blen, double dis_i, double dis_f);
