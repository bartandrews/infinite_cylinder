infinite_cylinder
=================

2D DMRG models on an infinite cylinder

Prerequisites: TeNPy, gnuplot

File naming convention
----------------------

All output .dat files are named in the following order:

*stem*

- tool (e.g. corr_len)
- charge (either "_charge" or nothing)
- model (e.g. Hubbard)
- lattice (e.g. Honeycomb)
- initial state (e.g. neel)
- tile units ([0, 1] or ['up', 'down'])
- chi (i.e. chi_max)

*leaf*

- t
- U
- mu
- V
- Lx
- Ly
- phi

NB: For a range of parameter values in an output file, we denote this by the order: min value _ max value _ number of samples (e.g. V_0_1_4).

Example: ent\_spec\_real\_charge\_Hubbard\_Square\_neel\_tile\_down\_up\_chi\_100\_t_\-1\_U\_1\_mu\_0.5\_V\_0\_Lx\_2\_Ly\_2.dat

Directory structure
-------------------

**data** is used to store all of the output dat files. The subdirectories are the output directories for the tools which I have defined (e.g. **ent_spec_real**). Inside each of the tools subdirectories there are the plotting scripts, as well as a **keep** subsubdirectory. It is intended that successful good-quality output is manually moved into keep. NB: No dat files are tracked by git due to their potentially large size.

**models** is used to store custom MPO Hamiltonian python class files. Basic Hamiltonians are already implemented in TeNPy (e.g. Ising model). However, in this directory we can create our own Hamiltonian classes e.g. for twisted bilayer graphene.

**.idea** is used to store PyCharm configuration files, in case I would like to make changes to the code using a PyCharm project on a remote computer.
