# N-body-analysis

This script computes clustering statistics (cross-correlation, auto-correlation, shell estimator and multipoles in Legendre basis) with various effects (gravitational redshift, light cone, special relativistic beaming, transverse Doppler effect, etc).

A slice plot for the simulation (starts from particle level):

<img src="https://cloud.githubusercontent.com/assets/6698757/19275240/15c923bc-8fa1-11e6-86d7-d93a75b8056f.png" width="500">

2D cross-correlation in redshift space:

<img src="https://cloud.githubusercontent.com/assets/6698757/19275218/01d4814e-8fa1-11e6-952c-0a6eb507435a.png" width="400">

## Getting started

```
$ git clone https://github.com/OMGitsHongyu/N-body-analysis.git
$ cd pkgs/; python setup.py build_ext --inplace
```

## Usage

The main scripts are `./script/nbody_mean.py` and `./script/nbody_shadab.py`. See bash script for submitting to the queue.
