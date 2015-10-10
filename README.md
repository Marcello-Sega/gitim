### Status
[![Build Status](https://travis-ci.org/Marcello-Sega/gitim.svg?branch=ITIM)](https://travis-ci.org/Marcello-Sega/gitim)

## What is ITIM/GITIM?

**ITIM** and **GITIM** are two algorithms for the identification
of interfacial molecules or atoms. GITIM is a generalization of the
ITIM algorithm, which was designed for planar interfaces. GITIM,
on the contrary, is free from any geometrical constraints and can
be used to analyze intrinsic properties of surfaces with arbitrary
shape. 

## What can I find here ? 

There are two branches, ITIM and GITIM, which implement the respective
algorithms in the form of a program that can read
[GROMACS](http://www.gromacs.org) trajectories.


## Which branch should I use, ITIM or GITIM? 

If you are interested in planar interfaces only, use the ITIM branch.

Currently the code is undergoing a major restructuring, with new
features added and bugs fixed. The ITIM branch is so far the most
up-to-date, with a new build system that ensures compatibility with
a broader range of GROMACS versions. New features/bugfix are being
backported to the GITIM branch.

## Which features are included? 


* Identification of interfacial atoms
* Identification of interfacial molecules
* Identification of further atomic/molecular layers
* Cluster search to deal with partially miscible systems
* Calculation of intrinsic and non-intrinsic density and mass profiles
* Use of enhanced Monte Carlo sampling to take into account the fluctuations of the opposite interface on the density profiles
* Exporting PDB files containing also the information on the layer to which atoms belong to, as well as the intrinsic distance from the surface.

_Note: Some of the features are so far present only in the ITIM branch and are currently being ported to the GITIM branch._


## Which code is faster, ITIM or GITIM

ITIM, by far. This is because GITIM needs to compute the alpha-shape
of the phase, which in turn relies on the Delauney triangulation.
This is done through QHull, but it is still much slower than the
ITIM algorithm.  So, if you only have to deal with planar interfaces,
go for the ITIM branch.

## What do I need to build the code ? 

* A gromacs build + sources (compatibility starting from GROMACS 4.6.7 up to 5.1) 
* Cmake >= 2.8.11

## What can you do with GITIM?

Check out the [usage example page](http://www.gitim.eu/usage-examples) of the gitim web site http://www.gitim.eu

## References

If you use this code to publish some research results, please read and cite both of the two following papers.

[M. Sega, S. S. Kantorovich P. Jedlovszky and M. Jorge, _J. Chem. Phys._ **138**, 044110 (2013)](http://dx.doi.org/10.1063/1.4776196) The generalized identification of truly interfacial molecules (ITIM) algorithm for nonplanar interfaces.

[L. B. Pártay, G. Hantal, P. Jedlovszky, Á. Vincze and G. Horvai, _J. Comp. Chem._ **29**, 945 (2008)] (http://dx.doi.org/10.1002/jcc.20852)
A new method for determining the interfacial molecules and characterizing the surface roughness in computer simulations. Application to the liquid–vapor interface of water

