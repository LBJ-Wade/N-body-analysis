#!/usr/bin/env python
"""
This program is used to compute monopole from 2d cross-correlation.
Usage: python file.py data gravfac tdfac rlcfac beamfac hxvfac xgfac masscut masserr wopec
"""

from __future__ import absolute_import, division
import os
import sys
import numpy as np
import time
import re
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(1, parent_dir)
from pkgs import correlation
from matplotlib import pylab as plt

__version__ = '2.0'
__author__ = 'Hongyu Zhu'

# -------------------------------- Define parameters -------------------------------------------------------
SPEED_OF_LIGHT = 299792.458
nbins = 25
nhocells = 50
blen = 1000
dis_i = 0.5
dis_f = 200
t0 = time.time()

# -------------------------------- Read from input ---------------------------------------------------------

print "Usage: python file.py data gravfac tdfac rlcfac beamfac hxvfac xgfac masscut masserr wopec"
gravfac = np.double(sys.argv[2])
tdfac = np.double(sys.argv[3])
rlcfac = np.double(sys.argv[4])
beamfac = np.double(sys.argv[5])
hxvfac = np.double(sys.argv[6])
xgfac = np.double(sys.argv[7])
masscut = np.double(sys.argv[8])
masserr = np.double(sys.argv[9])
wopec = np.double(sys.argv[10])

# -------------------------------- Output files ------------------------------------------------------------

outdir = '../est/newbeaming/test2/'

outfile_est = outdir + 'est_masscut_'+str(masscut)+'_masserr_'+str(masserr)+'_grav_'+sys.argv[2]+\
        '_td_'+sys.argv[3]+'_lc_'+sys.argv[4]+'_bm_'+sys.argv[5]+'_hxv_'+sys.argv[6]+'_xg_'+\
        sys.argv[7]+'_'+re.split('/|\.',sys.argv[1])[-2]+'.txt'
outfile_auto1 = outdir + 'auto1_masscut_'+str(masscut)+'_masserr_'+str(masserr)+'_grav_'+sys.argv[2]+\
        '_td_'+sys.argv[3]+'_lc_'+sys.argv[4]+'_bm_'+sys.argv[5]+'_hxv_'+sys.argv[6]+'_xg_'+\
        sys.argv[7]+'_'+re.split('/|\.',sys.argv[1])[-2]+'.txt'
outfile_auto2 = outdir + 'auto2_masscut_'+str(masscut)+'_masserr_'+str(masserr)+'_grav_'+sys.argv[2]+\
        '_td_'+sys.argv[3]+'_lc_'+sys.argv[4]+'_bm_'+sys.argv[5]+'_hxv_'+sys.argv[6]+'_xg_'+\
        sys.argv[7]+'_'+re.split('/|\.',sys.argv[1])[-2]+'.txt'
outfile_xc = outdir + 'xc_masscut_'+str(masscut)+'_masserr_'+str(masserr)+'_grav_'+sys.argv[2]+\
        '_td_'+sys.argv[3]+'_lc_'+sys.argv[4]+'_bm_'+sys.argv[5]+'_hxv_'+sys.argv[6]+'_xg_'+\
        sys.argv[7]+'_'+re.split('/|\.',sys.argv[1])[-2]+'.txt'
outfile_poles = outdir + 'poles_masscut_'+str(masscut)+'_masserr_'+str(masserr)+'_grav_'+sys.argv[2]+\
        '_td_'+sys.argv[3]+'_lc_'+sys.argv[4]+'_bm_'+sys.argv[5]+'_hxv_'+sys.argv[6]+'_xg_'+\
        sys.argv[7]+'_'+re.split('/|\.',sys.argv[1])[-2]+'.txt'

if wopec > 0:
    outfile_est = outfile_est[:-4]+'_pec.txt'
    outfile_auto1 = outfile_auto1[:-4]+'_pec.txt'
    outfile_auto2 = outfile_auto2[:-4]+'_pec.txt'
    outfile_xc = outfile_xc[:-4]+'_pec.txt'
    outfile_poles = outfile_poles[:-4]+'_pec.txt'

print "output estimator file:", outfile_est
print "output xc file:", outfile_xc
print "output auto1 file:", outfile_auto1
print "output auto2 file:", outfile_auto2
print "output xc monopole/dipole/quadrupole file:", outfile_poles

# ----------------------- Open the simulation data ----------------------------------------------------------

print "Mean potential is used"
galdata0 = np.loadtxt(sys.argv[1])
acc = np.loadtxt(sys.argv[11])
galdata = np.column_stack((galdata0[:,:8], -acc))
ngal = galdata.shape[0]

# ---------------------- Sort the particles by mass in a descending order ----------------------------------

sortedgal0 = galdata[np.array(-1.*galdata[:,7]).argsort()]
sortedgal = sortedgal0[sortedgal0[:,7] > masscut]
ngalcut = len(sortedgal)
print "Number of galaxies", ngal, "Number of galaxies left", ngalcut
sortedgal_mass = 10**(np.log10(sortedgal[:,7])+np.random.randn(sortedgal.shape[0])*masserr) # Add mass spread
sortedgal = sortedgal[np.array(-1.*sortedgal_mass).argsort()]

# ----------------------- Add effects ----------------------------------------------------------------------

probtokeepbm = 1.+(-sortedgal[:,5])*1*beamfac/SPEED_OF_LIGHT
weight = probtokeepbm

sortedgal = np.column_stack((sortedgal, weight))
sortedgal[:,0:3] = sortedgal[:,0:3] / 1000
sortedgal[:,6] = -sortedgal[:,6]/SPEED_OF_LIGHT
zz = sortedgal[:,2].copy()
vz = sortedgal[:,5].copy()

if gravfac > 0:
    print "gravitational redshift effect included with factor", gravfac
    sortedgal[:,2] += sortedgal[:,6] / 100 * gravfac
if tdfac > 0:
    print "transvers Doppler effect included with factor", tdfac
    for i in xrange(ngalcut):
        sortedgal[i,2] += (sortedgal[i,3]*sortedgal[i,3]+sortedgal[i,4]*sortedgal[i,4]\
                +sortedgal[i,5]*sortedgal[i,5]) / 100 * tdfac / SPEED_OF_LIGHT / 2
if rlcfac > 0:
    print "light cone effect included with factor", rlcfac
    for i in xrange(ngalcut):
        sortedgal[i,2] += (sortedgal[i,5]*sortedgal[i,5]) / 100 * rlcfac / SPEED_OF_LIGHT
if beamfac > 0:
    print "SRB effect included with factor", beamfac
if hxvfac > 0:
    print "Hxv effect included with factor", hxvfac
    sortedgal[:,2] += zz * vz / SPEED_OF_LIGHT * hxvfac
if xgfac > 0:
    print "xgfac effect included with factor", xgfac
    sortedgal[:,2] += - sortedgal[:,8] * zz * 1000 / SPEED_OF_LIGHT * xgfac / 100
#    sortedgal[:,2] += - sortedgal[:,8] * 1000 / SPEED_OF_LIGHT * xgfac / 100
if wopec > 0:
    sortedgal[:,2] += sortedgal[:,5] / 100
    print "peculiar velocities are included"

# ------------------------ Wrap the galaxies ---------------------------------------------------------------

for i in xrange(ngalcut):
    for j in xrange(3):
        sortedgal[i,j] = sortedgal[i,j] % blen
#        if sortedgal[i,j] < 0:
#            sortedgal[i,j] += blen
#        if sortedgal[i,j] >= blen:
#            sortedgal[i,j] -= blen

# ----------------------- Divide particles into 2 groups and make them c contiguous ------------------------
print np.where(np.cumsum(sortedgal[:,9]) >= sortedgal[:,9].sum()/2)
ngalcut2 = np.where(np.cumsum(sortedgal[:,9]) >= sortedgal[:,9].sum()/2)[0][0] + 1
g1 = np.ascontiguousarray(sortedgal[0:ngalcut2,[0,1,2,3,4,5,9]], dtype=np.float64)
g2 = np.ascontiguousarray(sortedgal[ngalcut2:,[0,1,2,3,4,5,9]], dtype=np.float64)
print "High mass sample:", g1.shape[0], g1[:,6].sum()
print "Low mass sample:", g2.shape[0], g2[:,6].sum()
print g1[0]

# ------------------------ Calculate the shell estimator ---------------------------------------------------

shell = correlation.estimator1dpy(g1, g2, nbins, nhocells, blen, dis_i, dis_f, 1)
np.savetxt(outfile_est, shell)

# ------------------------ Compute the correlation function ------------------------------------------------
"""
dis_f = 200
nbins = 40
rlim = 200
poles = correlation.polepy(g1, g2, rlim, nbins, nhocells, blen, dis_f)
np.savetxt(outfile_poles, poles)
poles1 = correlation.polepy(g1, g1, rlim, nbins, nhocells, blen, dis_f)
np.savetxt(outfile_auto1, poles1)
poles2 = correlation.polepy(g2, g2, rlim, nbins, nhocells, blen, dis_f)
np.savetxt(outfile_auto2, poles2)
nbins = 400
dis_f = 200
rlim = 50

dist0, dist1, xc, npair = correlation.corr2dpy(g1, g2, rlim, nbins, nhocells, blen, dis_f, 1, 0)
np.savetxt(outfile_xc, xc)
"""
print "\n", time.time() - t0
