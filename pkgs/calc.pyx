# f.pyx: numpy arrays -> extern from "fc.h"
# 3 steps:
# cython calc.pyx  -> f.c
# link: python f-setup.py build_ext --inplace  -> f.so, a dynamic library
# py test-f.py: import f gets f.so, f.fpy below calls fc()

from __future__ import division
import numpy as np
cimport numpy as np
from matplotlib import pyplot as plt

cdef extern from "./corr2dc.h" nogil:
    void corr2d(double *xc_2d, double *p1, long np1, double *p2, long np2, int nbins,
            int nhocells, double blen, double dis_f, int weighted, int ispolar)

    void estimator1d(double *xc_1d, double *npair, double *p1, long np1, double *p2, long np2,
            int nbins, int nhocells, double blen, double dis_i, double dis_f, int weighted);

    void estimator1dhaloterm(double *xc_1d, double *npair, double *p1, long np1, double *p2, long np2,
            int nbins, int nhocells, double blen, double dis_i, double dis_f, int weighted, int haloterm)

    void corr1d(double *xc_1d, double *p1, long np1, double *p2, long np2, int nbins,
            int nhocells, double blen, double dis_i, double dis_f);

def estimator1dpy(np.ndarray[np.double_t,ndim=2] p1,
        np.ndarray[np.double_t,ndim=2] p2, nbins, nhocells, blen, dis_i, dis_f,
        weighted):
    """
    Calculate the 1d shell estimator of 2d correlation function.

    Parameters
    ----------
    p1: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 1.
    p2: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 2.
    nbins: int
        Number of bins in the estimator.
    nhocells: int
        Number of cells.
    blen: double
        Box length.
    dis_i: double
        Minimum distance in the logspace.
    dis_f: double
        Maximum distance in the logspace.
    weighted: int
        For possible lightcone and beaming effects.
    Returns
    -------
    dist: (nbins,) array of double
        r of shell estimator.
    xc: (nbins,) array of double
        Shell estimator.
    npair: (nbins,) array of double
        Pair count.

    """
    cdef:
        np.ndarray[np.double_t,ndim=1] xc = np.zeros((nbins,), dtype=np.float64)
        np.ndarray[np.double_t,ndim=1] npair = np.zeros((nbins,), dtype=np.float64)
    nattr = 7 if weighted > 0 else 6
    assert p1.shape[1] == p2.shape[1] == nattr
    np1, np2 = p1.shape[0], p2.shape[0]
    rint = (np.log10(dis_f) - np.log10(dis_i)) / nbins
    dist = 10 ** (np.log10(dis_i) + (np.arange(nbins)+0.5)*rint)

    estimator1d(<double*> xc.data, <double*> npair.data, <double*> p1.data, np1, <double*> p2.data, np2,
            nbins, nhocells, blen, dis_i, dis_f, weighted)

    return np.column_stack((dist, xc, npair))


def estimator1d12py(np.ndarray[np.double_t,ndim=2] p1,
        np.ndarray[np.double_t,ndim=2] p2, nbins, nhocells, blen, dis_i, dis_f,
        weighted, nhaloterm):
    """
    Calculate the 1d shell estimator of 2d correlation function based on 1/2-halo term.

    Parameters
    ----------
    p1: (np1,7/8) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 1.
    p2: (np1,7/8) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 2.
    nbins: int
        Number of bins in the estimator.
    nhocells: int
        Number of cells.
    blen: double
        Box length.
    dis_i: double
        Minimum distance in the logspace.
    dis_f: double
        Maximum distance in the logspace.
    weighted: int
        For possible lightcone and beaming effects.
    nhaloterm: int
        n-halo term.
    Returns
    -------
    dist: (nbins,) array of double
        r of shell estimator.
    xc: (nbins,) array of double
        Shell estimator.
    npair: (nbins,) array of double
        Pair count.

    """
    cdef:
        np.ndarray[np.double_t,ndim=1] xc = np.zeros((nbins,), dtype=np.float64)
        np.ndarray[np.double_t,ndim=1] npair = np.zeros((nbins,), dtype=np.float64)
    nattr = 8 if weighted > 0 else 7
    assert p1.shape[1] == p2.shape[1] == nattr
    np1, np2 = p1.shape[0], p2.shape[0]
    rint = (np.log10(dis_f) - np.log10(dis_i)) / nbins
    dist = 10 ** (np.log10(dis_i) + (np.arange(nbins)+0.5)*rint)

    estimator1dhaloterm(<double*> xc.data, <double*> npair.data, <double*> p1.data, np1, <double*> p2.data, np2,
            nbins, nhocells, blen, dis_i, dis_f, weighted, nhaloterm)

    return np.column_stack((dist, xc, npair))

def corr1dpy(np.ndarray[np.double_t,ndim=2] p1,
        np.ndarray[np.double_t,ndim=2] p2, nbins, nhocells, blen, dis_i, dis_f):
    """
    Calculate the 1d correlation function. (e.g. dz vs r)

    Parameters
    ----------
    p1: (np1,4) array of double
        Particle x, y, z, feature of Group 1.
    p2: (np1,6/7) array of double
        Particle x, y, z, feature of Group 2.
    nbins: int
        Number of bins in the estimator.
    nhocells: int
        Number of cells.
    blen: double
        Box length.
    dis_i: double
        Minimum distance in the logspace.
    dis_f: double
        Maximum distance in the logspace.
    Returns
    -------
    xc: (nbins,) array of double
        1d correlation function.

    """
    cdef np.ndarray[np.double_t,ndim=1] xc = np.zeros((nbins,), dtype=np.float64)
    assert p1.shape[1] == p2.shape[1] == 4
    np1, np2 = p1.shape[0], p2.shape[0]

    corr1d(<double*> xc.data, <double*> p1.data, np1, <double*> p2.data, np2,
            nbins, nhocells, blen, dis_i, dis_f)

    return xc

def corr2dpy(np.ndarray[np.double_t,ndim=2] p1, np.ndarray[np.double_t,ndim=2] p2,
        rlim, nbins, nhocells, blen, dis_f, weighted, ispolar):
    """
    Calculate the 2d correlation function.

    Parameters
    ----------
    p1: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 1.
    p2: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 2.
    rlim: double
        Maximum distance in the xc/npair.
    nbins: int
        Number of bins in the estimator.
    nhocells: int
        Number of cells.
    blen: double
        Box length.
    dis_f: double
        Maximum distance.
    weighted: int
        For possible lightcone and beaming effects.
    ispolar : int
        Coordinate systems. Cartesian when ispolar <= 0.
    Returns
    -------
    xdist: (nlim-1,) array of double
        x of correlation function (r_perp, centered at 0).
        cos of correlation function if ispolar > 0.
    ydist: (nlim,) array of double
        y of correlation function (r_parallel).
        r of correlation function if ispolar > 0.
    xc: (nlim,nlim-1) array of double
        correlation function.
    npair: (nlim,nlim-1) array of double
        Pair count.

    """
    cdef np.ndarray[np.double_t,ndim=2] npair = np.zeros((nbins, nbins), dtype=np.float64)

    if weighted > 0:
        nattr = 7
        np1, np2 = p1[:,6].sum(), p2[:,6].sum()
    else:
        nattr = 6
        np1, np2 = p1.shape[0], p2.shape[0]

    assert p1.shape[1] == p2.shape[1] == nattr

    corr = np.zeros((nbins, nbins))
    r_perp = np.zeros((nbins+1,))
    simVolumn = blen ** 3

    corr2d(<double*> npair.data, <double*> p1.data, p1.shape[0], <double*> p2.data, p2.shape[0],
            nbins, nhocells, blen, dis_f, weighted, ispolar)

    if ispolar > 0:

        rint = dis_f / nbins
        nlim, n_f = np.floor(rlim / rint), np.floor(dis_f / rint)
        for i in xrange(len(r_perp)):
            r_perp[i] = rint * i
        for i in xrange(nbins):
            corr[:,i] = 0.75 * npair[:,i] / np1 / np2 * simVolumn / np.pi /\
           (r_perp[i+1]**3 - r_perp[i]**3) * nbins - 1

        xdist = (np.arange(1,nlim) - 0.5) / nbins
        ydist = (np.arange(nlim) + 0.5) * rint
        xc = corr[:,0:nlim]
        xc_npair = npair[:,0:nlim]

    else:

        rint = 2 * dis_f / nbins
        nlim, n_f = np.floor(rlim / rint), np.floor(dis_f / rint)
        for i in xrange(len(r_perp)):
            r_perp[i] = rint * i
        for i in xrange(nbins):
            corr[i,:] = npair[i,:] / np1 / np2 * simVolumn / blen / np.pi /\
           (r_perp[i+1]**2 - r_perp[i]**2) * nbins - 1

        xdist = (np.arange(1,2*nlim) - nlim) * rint
        ydist = (np.arange(0,2*nlim) - nlim + 0.5) * rint
        xc = np.transpose(np.row_stack((corr[::-1][:nbins-1,:], corr)))[n_f-nlim:n_f+nlim, nbins-nlim:nbins+nlim-1]
        xc_npair = np.transpose(np.row_stack((npair[::-1][:nbins-1,:], npair)))[n_f-nlim:n_f+nlim, nbins-nlim:nbins+nlim-1]

    return xdist, ydist, xc, xc_npair


def monopolepy(np.ndarray[np.double_t,ndim=2] p1,
        np.ndarray[np.double_t,ndim=2] p2, rlim, nbins, nhocells, blen, dis_f):
    """
    Calculate the monopole of 2d correlation function.

    Parameters
    ----------
    p1: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 1.
    p2: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 2.
    rlim: double
        Maximum range to look at.
    nbins: int
        Number of bins in the estimator.
    nhocells: int
        Number of cells.
    blen: double
        Box length.
    dis_f: double
        Maximum distance.
    Returns
    -------
    monopole_obs: (nbins,2) array of long
        Monopole of 2d correlation function.

    """
    cdef np.ndarray[np.double_t,ndim=2] npair = np.zeros((nbins, nbins), dtype=np.float64)

    assert p1.shape[1] == p2.shape[1]
    np1, np2 = p1.shape[0], p2.shape[0]

    corr = np.zeros((nbins, nbins))
    r_perp = np.zeros((nbins+1,))
    r = np.zeros((nbins,))
    simVolumn = blen ** 3

    corr2d(<double*> npair.data, <double*> p1.data, np1, <double*> p2.data, np2,
            nbins, nhocells, blen, dis_f, 1, 1)

    rint = dis_f / nbins
    for i in xrange(len(r_perp)):
        r_perp[i] = rint * i
    for i in xrange(nbins):
        corr[:,i] = 0.75 * npair[:,i] / np1 / np2 * simVolumn / np.pi /\
       (r_perp[i+1]**3 - r_perp[i]**3) * nbins - 1

    monopole = np.zeros((nbins,))
    for i in xrange(nbins):
        r[i] = rint * (i+0.5)
        for j in xrange(nbins):
            monopole[i] += corr[j,i] * 2. / nbins

    monopole_obs = np.column_stack((r, monopole))

    return monopole_obs


def polepy(np.ndarray[np.double_t,ndim=2] p1,
        np.ndarray[np.double_t,ndim=2] p2, rlim, nbins, nhocells, blen, dis_f):
    """
    Calculate the quadrupole of 2d correlation function.

    Parameters
    ----------
    p1: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 1.
    p2: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 2.
    rlim: double
        Maximum range to look at.
    nbins: int
        Number of bins in the estimator.
    nhocells: int
        Number of cells.
    blen: double
        Box length.
    dis_f: double
        Maximum distance.
    Returns
    -------
    pole_obs: (nbins,4) array of double
        r, monopole, dipole, quadrupole of 2d correlation function.

    """
    cdef np.ndarray[np.double_t,ndim=2] npair = np.zeros((nbins, nbins), dtype=np.float64)

    assert p1.shape[1] == p2.shape[1]
#    np1, np2 = p1.shape[0], p2.shape[0]
    np1, np2 = p1[:,6].sum(), p2[:,6].sum()

    corr = np.zeros((nbins, nbins))
    r_perp = np.zeros((nbins+1,))
    r, mu = np.zeros((nbins,)), np.zeros((nbins,))
    simVolumn = blen ** 3

    corr2d(<double*> npair.data, <double*> p1.data, p1.shape[0], <double*> p2.data, p2.shape[0],
            nbins, nhocells, blen, dis_f, 1, 1)

    rint = dis_f / nbins
    for i in xrange(len(r_perp)):
        r_perp[i] = rint * i
    for i in xrange(nbins):
        corr[:,i] = 0.75 * npair[:,i] / np1 / np2 * simVolumn / np.pi /\
       (r_perp[i+1]**3 - r_perp[i]**3) * nbins - 1

    monopole, dipole, quadrupole = np.zeros((nbins,)), np.zeros((nbins,)), np.zeros((nbins,))

    mu = np.linspace(-1, 1, nbins+1)
    dmu = 2 / nbins
    for i in xrange(nbins):
        r[i] = rint * (i+0.5)
        for j in xrange(nbins):
            monopole[i] += corr[j,i] * dmu / 2
            dipole[i] += corr[j,i] * mu[j] * dmu * 3 / 2
            quadrupole[i] += corr[j,i] * (3*mu[j]*mu[j]-1) * dmu / 4 * 5

    pole_obs = np.column_stack((r, monopole, dipole, quadrupole))

    return pole_obs
