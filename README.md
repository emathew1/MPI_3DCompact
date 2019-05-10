# MPI_3DCompact

CFD framework for compressible flow solvers using compact finite differences and Runge-Kutta timestepping methods. Parallelized using MPI and 2-D pencil domain decomposition.  

Current default (and only up to date) solver is the `CurvilinearCSolver`, solves compressible Navier-Stokes equations on a single, generalized curvlinear block using the formulation of Visbal & Gaitonde (2002, JCP). There's also a uniform grid, compressible solver with fully-expanded viscous terms that needs to be brought up to speed in `DepreciatedSolvers.` 

## Required Libraries
*libpng* - required for generating flow-field images

*boost/program_options* - required for parsing flow solver options inputs

*whatever mpi you'd like* - if you like running things in parallel

## Spatial Schemes 

`*FDTYPES` = `CD2, PADE6, PENTA10`


`CD2` - Second-order Central Difference Scheme

`PADE6` - Sixth-order Tridiagonal Compact Pade Scheme (Lele, 1992)

`PENTA10` - Tenth-order Pentadiagonal Compact Scheme (Lele, 1992)

## Temporal Schemes
 `RKTYPE` = `TVDRK3, RK4, KENRK4, LSLDDRK4`
 
 
 `TVDRK3` - Total-Variation Diminishing RK3 scheme (Shu & Osher, 1988)
 
 `RK4` - Classic Fourth-order RK method 
 
 `KENRK4` - Five-stage Fourth-order Lowish-Dispersion RK4 (Kennedy et al, 1999)
 
 `LSLDDRK4` - Two Stage, Eleven Step Low-Storage, Low-Dispersion-Dissipation, RK4 (Stanescu & Habashi, 1998)
 
## Filters

`FILTERTYPE` = `COMPACT8, COMPACT10`

`COMPACT8` - Eigth-order compact filter (Giatonde & Visbal, 2000)

`COMPACT10` - Tenth-order compact filter (Giatonde & Visbal, 2000)

## LES Models

`LESMODEL` = `NONE, VREMAN`

`NONE` = The special DNS LES model

`VREMAN` = Standard Vreman SGS (Vreman, 2004)

`LESAVERAGING` = `NONE, XI1_AVG, XI2_AVG, XI3_AVG, LOCAL`

`XI1_AVG, XI2_AVG, XI3_AVG` - Averaging along cartesian (pre-transformed) directions

`LOCAL` - Using a filter function to do a local filtering/averaging of quantities

