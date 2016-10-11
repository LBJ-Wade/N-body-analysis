#include "corr2dc.h"

void init_mesh(long *ll, long *hoc, double *p, long np, int nattr, int nho, double blen) {

    // Initialize the space. Mesh the particles in p into cells.
    //
    // Parameters:
    // ll : (np,) array, linked list of particles.
    // hoc : (nhocells, nhocells, nhocells) array, mesh of particles.
    // p : (np, nattr) array, particles.
    // np : number of particles.
    // nattr: number of features of particles.
    // nho : number of cells in one dimension.
    // blen : box length.

    long i;
    int ix, iy, iz, ii, jj, kk;

    // Initiate hoc to -1 for not missing the last point
    for (ii = 0; ii < nho; ii++)
        for (jj = 0; jj < nho; jj++)
            for (kk = 0; kk < nho; kk++)
                hoc[ii*nho*nho+jj*nho+kk] = -1;

    for (i = 0; i < np; i++) {
        ix = floor(p[i*nattr] / blen * nho);
        iy = floor(p[i*nattr+1] / blen * nho);
        iz = floor(p[i*nattr+2] / blen * nho);

        ll[i] = hoc[ix*nho*nho+iy*nho+iz];
        hoc[ix*nho*nho+iy*nho+iz] = i;
    }
}

void corr2d(double *xc_2d, double *p1, long np1, double *p2, long np2, int nbins,
        int nhocells, double blen, double dis_f, int weighted, int ispolar) {

    // 2d correlation function.
    //
    // Parameters:
    // xc_2d : (nbins, nbins) array, 2d correlation to be returned.
    // p1 : (np1,6/7) first group of particles.
    // np1: number of particles in first group.
    // p2 : (np2,6/7) second group of particles.
    // np2 : number of particles of the second group.
    // nbins: number of bins.
    // nhocells: number of cells.
    // dis_f : the maximum range in xc_2d.
    // weighted : this is due to the possible lightcone and beaming effect.
    // ispolar : coordinate system, polar/cartesian.

    long iq, i;
    int pp, qq, rr, p, q, r, ix, iy, iz, k_perp, k_para, k_r, k_costheta;
    double xqi, yqi, zqi, dx, dy, dz, dr, wt;
    int nattr = 6, ncb = floor((dis_f / blen) * (double)(nhocells)) + 1;

    long *ll = (long *)calloc(np2, sizeof(long));
    long *hoc = (long *)calloc(nhocells*nhocells*nhocells, sizeof(long));

    // Add weight or not
    if (weighted > 0)  nattr = 7;

    // initialization
    init_mesh(ll, hoc, p2, np2, nattr, nhocells, blen);

    for (iq = 0; iq < np1; iq++) {
        if (iq % 500 == 0) printf("%ld", iq);

        xqi = p1[iq*nattr];
        yqi = p1[iq*nattr+1];
        zqi = p1[iq*nattr+2];

        ix = floor(xqi / blen * nhocells);
        iy = floor(yqi / blen * nhocells);
        iz = floor(zqi / blen * nhocells);

        for (pp = ix - ncb; pp < ix + ncb; pp++) {
            p = CWRAP(pp, nhocells);
            for (qq = iy - ncb; qq < iy + ncb; qq++) {
                q = CWRAP(qq, nhocells);
                for (rr = iz - ncb; rr < iz + ncb; rr++) {
                    r = CWRAP(rr, nhocells);
                    if (hoc[p*nhocells*nhocells+q*nhocells+r] != -1) {
                        i = hoc[p*nhocells*nhocells+q*nhocells+r];
                        while (1) {

                            dx = RELWRAP(p2[i*nattr] - xqi, blen);
                            dy = RELWRAP(p2[i*nattr+1]- yqi, blen);
                            dz = RELWRAP(p2[i*nattr+2]- zqi, blen);

                            wt = 1;
                            if (weighted > 0)  wt = p1[iq*nattr+6] * p2[i*nattr+6];

                            if (ispolar > 0) {

                                dr = sqrt(dx * dx + dy * dy + dz * dz);
                                k_r = floor(dr / dis_f * nbins);
//                                k_costheta = floor(sqrt(dx * dx + dy * dy) / dr * nbins);
                                if (dz < 0) {
                                    k_costheta = floor((1 - sqrt(dz * dz) / dr) * nbins / 2);
                                } else {
                                    k_costheta = floor((1 + sqrt(dz * dz) / dr) * nbins / 2);
                                }

                                // if there is any outlier in k_r
                                // if there is any outlier in k_costheta
                                if (k_r >= 0 && k_costheta >= 0 && k_r < nbins && k_costheta < nbins)
                                    xc_2d[k_costheta*nbins+k_r] += wt;

                            } else {

                                k_para = floor((dz + dis_f) * nbins / 2 / dis_f);
                                // if there is any outlier in k_para
                                if (k_para >= 0 && k_para < nbins) {
                                    dr = sqrt(dx * dx + dy * dy);
                                    k_perp = floor(dr * nbins / 2 / dis_f);
                                    xc_2d[k_perp*nbins+k_para] += wt;
                                }
                            }

                            if (ll[i] != -1) {
                                i = ll[i];
                            }
                            else break;
                        }
                    }
                }
            }
        }
    }
    free(ll);
    free(hoc);
}

void estimator1d(double *xc_1d, double *npair, double *p1, long np1, double *p2, long np2,
        int nbins, int nhocells, double blen, double dis_i, double dis_f, int weighted) {

    // 1d estimator of 2d correlation function.
    //
    // Parameters:
    // xc_1d : (nbins,) array, 1d estimator of 2d correlation to be returned.
    // p1 : (np1,6/7) first group of particles.
    // np1: number of particles in first group.
    // p2 : (np2,6/7) second group of particles.
    // np2 : number of particles of the second group.
    // nbins: number of bins.
    // nhocells: number of cells.
    // dis_i : since it's in the log space, specify the minimum range in xc_1d.
    // dis_f : the maximum range in xc_1d.
    // weighted : this is due to the possible lightcone and beaming effect.

    long iq, i;
    int pp, qq, rr, p, q, r, ix, iy, iz, k;
    double xqi, yqi, zqi, dx, dy, dz, dr, wt;
    int nattr = 6, ncb = floor((dis_f / blen) * (double)(nhocells)) + 1;
    double logdis_i = log10(dis_i);
    double logdis_f = log10(dis_f);

    // Add weight or not
    if (weighted > 0)  nattr = 7;

    long *ll = (long *)calloc(np2, sizeof(long));
    long *hoc = (long *)calloc(nhocells*nhocells*nhocells, sizeof(long));
    double *sumdz = (double *)calloc(nbins, sizeof(double));

    init_mesh(ll, hoc, p2, np2, nattr, nhocells, blen);

    for (iq = 0; iq < np1; iq++) {
        if (iq % 500 == 0) printf("%ld", iq);

        xqi = p1[iq*nattr];
        yqi = p1[iq*nattr+1];
        zqi = p1[iq*nattr+2];

        ix = floor(xqi / blen * nhocells);
        iy = floor(yqi / blen * nhocells);
        iz = floor(zqi / blen * nhocells);

        for (pp = ix - ncb; pp < ix + ncb; pp++) {
            p = CWRAP(pp, nhocells);
            for (qq = iy - ncb; qq < iy + ncb; qq++) {
                q = CWRAP(qq, nhocells);
                for (rr = iz - ncb; rr < iz + ncb; rr++) {
                    r = CWRAP(rr, nhocells);
                    if (hoc[p*nhocells*nhocells+q*nhocells+r] != -1) {
                        i = hoc[p*nhocells*nhocells+q*nhocells+r];
                        while (1) {

                            dx = RELWRAP(p2[i*nattr] - xqi, blen);
                            dy = RELWRAP(p2[i*nattr+1]- yqi, blen);
                            dz = RELWRAP(p2[i*nattr+2]- zqi, blen);

                            wt = 1;
                            if (weighted > 0)  wt = p1[iq*nattr+6] * p2[i*nattr+6];

                            dr = sqrt(dx * dx + dy * dy + dz * dz);
                            k = floor((log10(dr) - logdis_i) * nbins / (logdis_f - logdis_i));

                            // outliers
                            if (k >= 0 && k < nbins) {
                                npair[k] += wt;
                                sumdz[k] += dz * wt;
                            }

                            if (ll[i] != -1) {
                                i = ll[i];
                            }
                            else break;
                        }
                    }
                }
            }
        }
    }
    for (k = 0; k < nbins; k++) {
        if (npair[k] > 0) {
            xc_1d[k] = sumdz[k] / npair[k];
        } else {
            xc_1d[k] = 0;
        }
    }
    free(sumdz);
    free(ll);
    free(hoc);
}


void estimator1dhaloterm(double *xc_1d, double *npair, double *p1, long np1, double *p2, long np2,
        int nbins, int nhocells, double blen, double dis_i, double dis_f, int weighted, int haloterm) {

    // 1d estimator of 2d correlation function.
    //
    // Parameters:
    // xc_1d : (nbins,) array, 1d estimator of 2d correlation to be returned.
    // p1 : (np1,7/8) first group of particles.
    // np1: number of particles in first group.
    // p2 : (np2,7/8) second group of particles.
    // np2 : number of particles of the second group.
    // nbins: number of bins.
    // nhocells: number of cells.
    // dis_i : since it's in the log space, specify the minimum range in xc_1d.
    // dis_f : the maximum range in xc_1d.
    // weighted : this is due to the possible lightcone and beaming effect.

    long iq, i;
    int pp, qq, rr, p, q, r, ix, iy, iz, k, flag;
    double xqi, yqi, zqi, bqi, dx, dy, dz, dr, wt;
    int nattr = 7, ncb = floor((dis_f / blen) * (double)(nhocells)) + 1;
    double logdis_i = log10(dis_i);
    double logdis_f = log10(dis_f);

    // Add weight or not
    if (weighted > 0)  nattr = 8;

    long *ll = (long *)calloc(np2, sizeof(long));
    long *hoc = (long *)calloc(nhocells*nhocells*nhocells, sizeof(long));
    double *sumdz = (double *)calloc(nbins, sizeof(double));

    init_mesh(ll, hoc, p2, np2, nattr, nhocells, blen);

    for (iq = 0; iq < np1; iq++) {
        if (iq % 500 == 0) printf("%ld", iq);

        xqi = p1[iq*nattr];
        yqi = p1[iq*nattr+1];
        zqi = p1[iq*nattr+2];
        bqi = p1[iq*nattr+7];

        ix = floor(xqi / blen * nhocells);
        iy = floor(yqi / blen * nhocells);
        iz = floor(zqi / blen * nhocells);

        for (pp = ix - ncb; pp < ix + ncb; pp++) {
            p = CWRAP(pp, nhocells);
            for (qq = iy - ncb; qq < iy + ncb; qq++) {
                q = CWRAP(qq, nhocells);
                for (rr = iz - ncb; rr < iz + ncb; rr++) {
                    r = CWRAP(rr, nhocells);
                    if (hoc[p*nhocells*nhocells+q*nhocells+r] != -1) {
                        i = hoc[p*nhocells*nhocells+q*nhocells+r];

                        while (1) {

                            if (haloterm == 1) {
                                flag = (p2[i*nattr+7] == bqi);
                            } else if (haloterm == 2) {
                                flag = (p2[i*nattr+7] != bqi);
                            } else {
                                flag = 1;
                            }

                            if (flag) {

                                dx = RELWRAP(p2[i*nattr] - xqi, blen);
                                dy = RELWRAP(p2[i*nattr+1]- yqi, blen);
                                dz = RELWRAP(p2[i*nattr+2]- zqi, blen);

                                wt = 1;
                                if (weighted > 0)  wt = p1[iq*nattr+6] * p2[i*nattr+6];

                                dr = sqrt(dx * dx + dy * dy + dz * dz);
                                k = floor((log10(dr) - logdis_i) * nbins / (logdis_f - logdis_i));

                                // outliers
                                if (k >= 0 && k < nbins) {
                                    npair[k] += wt;
                                    sumdz[k] += dz * wt;
                                }
                            }

                            if (ll[i] != -1) {
                                i = ll[i];
                            }
                            else break;
                        }
                    }
                }
            }
        }
    }
    for (k = 0; k < nbins; k++) {
        if (npair[k] > 0) {
            xc_1d[k] = sumdz[k] / npair[k];
        } else {
            xc_1d[k] = 0;
        }
    }
    free(sumdz);
    free(ll);
    free(hoc);
}

void corr1d(double *xc_1d, double *p1, long np1, double *p2, long np2,
        int nbins, int nhocells, double blen, double dis_i, double dis_f) {

    // 1d correlation function. (e.g. z_g vs. r)
    //
    // Parameters:
    // xc_1d : (nbins,) array, 1d correlation to be returned.
    // p1 : (np1,6/7) first group of particles.
    // np1: number of particles in first group.
    // p2 : (np2,6/7) second group of particles.
    // np2 : number of particles of the second group.
    // nbins: number of bins.
    // nhocells: number of cells.
    // dis_i : since it's in the log space, specify the minimum range in xc_1d.
    // dis_f : the maximum range in xc_1d.
    // weighted : this is due to the possible lightcone and beaming effect.

    long iq, i;
    int pp, qq, rr, p, q, r, ix, iy, iz, k;
    double xqi, yqi, zqi, vqi, dx, dy, dz, dv, dr;
    int nattr = 4, ncb = floor((dis_f / blen) * (double)(nhocells)) + 1;
    double logdis_i = log10(dis_i);
    double logdis_f = log10(dis_f);


    long *ll = (long *)calloc(np2, sizeof(long));
    long *hoc = (long *)calloc(nhocells*nhocells*nhocells, sizeof(long));

    long *npair = (long *)calloc(nbins, sizeof(long));
    double *sumdz = (double *)calloc(nbins, sizeof(double));

    init_mesh(ll, hoc, p2, np2, nattr, nhocells, blen);

    for (iq = 0; iq < np1; iq++) {
        if (iq % 500 == 0) printf("%ld", iq);

        xqi = p1[iq*nattr];
        yqi = p1[iq*nattr+1];
        zqi = p1[iq*nattr+2];
        vqi = p1[iq*nattr+3];

        ix = floor(xqi / blen * nhocells);
        iy = floor(yqi / blen * nhocells);
        iz = floor(zqi / blen * nhocells);

        for (pp = ix - ncb; pp < ix + ncb; pp++) {
            p = CWRAP(pp, nhocells);
            for (qq = iy - ncb; qq < iy + ncb; qq++) {
                q = CWRAP(qq, nhocells);
                for (rr = iz - ncb; rr < iz + ncb; rr++) {
                    r = CWRAP(rr, nhocells);
                    if (hoc[p*nhocells*nhocells+q*nhocells+r] != -1) {
                        i = hoc[p*nhocells*nhocells+q*nhocells+r];
                        while (1) {

                            dx = RELWRAP(p2[i*nattr] - xqi, blen);
                            dy = RELWRAP(p2[i*nattr+1]- yqi, blen);
                            dz = RELWRAP(p2[i*nattr+2]- zqi, blen);
                            dv = p2[iq*nattr+3] - vqi;

                            dr = sqrt(dx * dx + dy * dy + dz * dz);
                            k = floor((log10(dr) - logdis_i) * nbins / (logdis_f - logdis_i));

                            // outliers
                            if (k >= 0 && k < nbins) {
                                npair[k] += 1;
                                sumdz[k] += dv;
                            }

                            if (ll[i] != -1) {
                                i = ll[i];
                            }
                            else break;
                        }
                    }
                }
            }
        }
    }
    for (k = 0; k < nbins; k++) {
        if (npair[k] > 0) {
            xc_1d[k] = sumdz[k] / npair[k];
        } else {
            xc_1d[k] = 0;
        }
    }
    free(npair);
    free(sumdz);
    free(ll);
    free(hoc);
}
