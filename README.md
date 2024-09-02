# SUNSET code

The **S**calable, **U**nstructured **N**ode-**SET** code, for DNS of compressible reacting flows.

Developed by Dr Jack King, University of Manchester.

## Overview and features

- Numerically solves the compressible Navier-Stokes equations.
- Domain discretised with unstructured node-set.
- Spatial discretisation between 4th and 10th order using **LABFM**.
- Temporal discretisation 3rd order explicit Runge-Kutta.
- Characteristic based boundary conditions
   + Walls (can be curved)
   + Inflow (subsonic, non-reflecting or hard)
   + Outflow (subsonic, non-reflecting)
   + Symmetry
   + Periodic
- Parallelised with OpenMP + MPI.
- Three-dimensional domain discretisation, tested up to 1024 cores.

## References

- LABFM fundamentals: King, Lind & Nasar (2020) JCP 415:109549 https://doi.org/10.1016/j.jcp.2020.109549
- LABFM fundamentals + BCs: King & Lind (2022) JCP 449:110760 https://doi.org/10.1016/j.jcp.2021.110760
- LABFM/FD for combustion: King (2023) https://arxiv.org/abs/2310.02200


