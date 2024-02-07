# sunset-flames
High-order mesh-free DNS for combustion simulations

The **S**calable **U**nstructured **N**ode-**SET** code for combustion

This code solves the fully compressible multi-species Navier-Stokes equations. 

Features:
  - Spatial discretisation with LABFM, between 4th and 10th order (default 8).
  - Temporal integration with RK3 and PID controlled time-stepping.
  - Inflow, wall and outflow boundaries with characteristic boundary condition formulation.
  - Periodic, symmetric, or 2nd order walls also available at external domain boundaries.
  - Spatially varying resolution

References:
  - King, Lind & Nasar (2020) JCP
  - King & Lind (2022) JCP
  - King (2024) CMAME

