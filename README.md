# FDFD_modules
personal set of FDFD codes used for my research. Much of this has been inspired by Yu (Jerry) Shi and Wonseok Shin

## dependencies
Functional version of MatLab 2016 or greater. This was coded in Matlab 2017.

## dimensionalities d = 1,2,3
Implementations of FDFD for 1D, 2D, and 3D simulations

## driven simulations and band structures
Will do simulations with a source (point source, line source, mode profiles, plane waves, etc.) or can do a variety of eigenmode solves for band structures

## Dispersive 2D mode-solvers/eigen solvers
New type of mode solvers which can be used for 2 dimensional structures, based on the paper shown here
https://www.osapublishing.org/ol/abstract.cfm?uri=ol-40-6-1053

As an example, we consider the classic and simple 2D photonic crystal circle (pillars) in the TE polarization (the parameters of the circle and unit cell are taken from the Johannopoulos book, Photonic Crystals).

![](img/TE_benchmarking_PWEM_and_FDFD_dispersive.png?raw=true)
Here, blue is FDFD and green is PWEM (need to fix the legend). There is an anomalous point which could have been caught with a mode filter (note no mode filter applied so almost everything output is realistic)

For comparison, here is the non-dispersive eigensolver (k is the input, frequency is the output)

![](img/TE_benchmarking_PWEM_and_FDFD_nondispsersive.png?raw=true)

The advantage of this formulation is you get the imaginary parts of the band structure 
## 3D Solver with acceleration and enhancements
uses a powerful reformulation of Maxwell's equations to accelerate iterative solutions based on the Beltrami-Laplace operator shown here: http://www.mit.edu/~wsshin/pdf/shin2013oe.pdf

tested on 60x60x60 grids (dipole in vacuum) on a laptop and QMR and runs reasonably quick

## basic adjoints
Very simple adjoint example (hopefully more work here later)

### unit conventions
For this package, note that most things are specified so that any spatial units are in units of microns.
