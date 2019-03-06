infinite_cylinder
=================

2D DMRG models on an infinite cylinder

Prerequisites: TeNPy, gnuplot

File naming convention
----------------------

All output .dat files are named in the following order:

- tool (e.g. corr_len)
- model (e.g. Hubbard)
- lattice (e.g. Honeycomb)
- initial state (e.g. neel)
- tile units (0 1 or up down)
- t
- U
- mu
- V
- Lx
- Ly
- phi

NB: For a range of parameter values in an output file, we denote this by the order: min value _ max value _ number of samples (e.g. V_0_1_4)
