infinite_cylinder
=================

2D DMRG models on an infinite cylinder

Prerequisites: TeNPy 0.3+, gnuplot, python 3+

Workflow
--------

The code is split into three independent parts to optimize performance.

**phi_flow** is a program that smoothly varies the external flux through the cylinder. This is used to identify a topological phase and calculate the Chern number.

**Ly_flow** is a program that varies the circumference length of the cylinder. This is used to calculate topological entropy and CFT edge counting.

**V_flow** is a program that varies the off-site interaction strength, as defined in the model Hamiltonian. This is used to characterize phase transitions e.g. metal to FCI.

The tools employed for each 'flow' are given in the table below.

========   ========
**flow**   **tool**
========   ========
phi_flow   * charge_pump
           * ent_spec_flow
--------   --------
Ly_flow    * ent_scal
           * ent_spec_real
           * ent_spec_mom
--------   --------
V_flow     * corr_len
           * ent_spec_V_flow
========   ========

Tools description
-----------------

The initial tool set is inspired by the paper: "Characterization and stability of a fermionic ν=1/3 fractional Chern insulator" by Adolfo G. Grushin, Johannes Motruk, Michael P. Zaletel, Frank Pollmann, PRB **91**, 035136 (2015). https://arxiv.org/abs/1407.6985

* charge_pump = charge pump

    This function is designed to plot the equivalent of Figs. 2.a,c) from the paper.

* ent_spec_flow = entanglement spectrum flow

    This function is designed to plot the equivalent of Figs. 2.b,d) from the paper.

* ent_scal = entanglement scaling

    This function is designed to plot the equivalent of Fig. 3.a) from the paper.

* ent_spec_real = entanglement spectrum in real space

    Entanglement energy as a function of bond in the unit cell.

* ent_spec_mom = entanglement spectrum in momentum space

    This function is designed to plot the equivalent of Fig. 3.b) from the paper.

* corr_len = correlation length

    This function is designed to plot the equivalent of the inset in Fig. 3.c) from the paper.

* ent_spec_V_flow = entanglement spectrum flow with respect to V

    This function is designed to plot the equivalent of Fig. 3.c) from the paper.

Models description
------------------

* bosons_haldane = hardcore boson Haldane model

    Equation 1 of "Characterizing topological order by studying the ground states of an infinite cylinder" by Lukasz Cincio and Guifre Vidal, PRL **110**, 067208 (2013). https://arxiv.org/abs/1208.2623

* fermions_haldane = spinless fermion Haldane model

    Equation 1 of "Characterization and stability of a fermionic ν=1/3 fractional Chern insulator" by Adolfo G. Grushin, Johannes Motruk, Michael P. Zaletel, Frank Pollmann, PRB **91**, 035136 (2015). https://arxiv.org/abs/1407.6985

* fermions_TBG1 = spinful fermions with two orbitals, Hubbard model for twisted bilayer graphene

    Section IV of "Model for the metal-insulator transition in graphene superlattices and beyond" by Noah F. Q. Yuan and Liang Fu, PRB **98**, 045103 (2018). https://arxiv.org/abs/1803.09699

* fermions_TBG2 = spinless fermions with two orbitals, tight-binding model for twisted bilayer graphene

    Section III of "Model for the metal-insulator transition in graphene superlattices and beyond" by Noah F. Q. Yuan and Liang Fu, PRB **98**, 045103 (2018). https://arxiv.org/abs/1803.09699



Directory structure
-------------------

**data** is used to store all of the output dat files, organised into their corresponding subdirectories. The subdirectories are the output directories for the tools which I have defined (e.g. **ent_spec_real**). Inside each of the tools subdirectories there are the plotting scripts, as well as a **keep** subsubdirectory. It is intended that successful good-quality output is manually moved into keep. NB: No dat files are tracked by git due to their potentially large size.

**code** contains the source code, split into the three independent parts: phi_flow, Ly_flow, and V_flow. **code/models** is used to store custom MPO Hamiltonian python class files. Basic Hamiltonians are already implemented in TeNPy (e.g. Ising model). However, in this directory we store our own Hamiltonian classes e.g. for twisted bilayer graphene.

**scripts** contains all of the SLURM batch scripts used for Hydra and Piz Daint.

**logs** is used to store all of the stdout and stderr files from the Hydra and Piz Daint batch scripts. NB: No log files are tracked by git.

**.idea** is used to store PyCharm configuration files, in case I would like to make changes to the code using a PyCharm project on a remote computer.

File naming convention
----------------------

All output .dat files are named in the following order:

*stem*

- tool (e.g. corr_len)
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
- phi (i.e. phi_ext)

NB: For a range of parameter values in an output file, we denote this by the order: min value _ max value _ number of samples (e.g. V_0_1_4). Or for discrete parameters, simply by the order: min value _ max value (e.g. Ly_2_8).

*name = stem + leaf*

Example:  data/ent\_spec\_real/ent\_spec\_real\_Hubbard\_Square\_neel\_tile\_down\_up\_chi\_100\_t_\-1\_U\_1\_mu\_0.5\_V\_0\_Lx\_2\_Ly\_2.dat

Benchmarks
----------

These benchmarks will be updated as the code is optimized. The specifications of the machines is as follows:

**laptop** is a Dell XPS13 9370 with a quad-core Intel Core i7-8550U CPU @ 1.80GHz and 16GB of RAM. One node with 4 OpenMP threads.

**Hydra** specification is listed here: https://www.id.uzh.ch/en/scienceit/infrastructure/hydra.html Flows spread across 3 nodes with 16 OpenMP threads each.

**Piz Daint** specification is listed here: https://www.cscs.ch/computers/piz-daint/ Flows spread across 3 nodes with 36 OpenMP threads each.

NB: TeNPy uses Intel MKL and so is OpenMP capable. However, neither MPI nor hyperthreading are configured at the moment.

----

1) half-filled Haldane Chern insulator (chi = 100, Ly = 3)

=========   ========   ========   ========
machine     phi_flow   Ly_flow    V_flow
=========   ========   ========   ========
laptop      00:00:00   00:02:08   00:00:00
---------   --------   --------   --------
Hydra       00:32:27   00:05:23   00:18:38
---------   --------   --------   --------
Piz Daint   00:00:00   00:00:00   00:00:00
=========   ========   ========   ========

2) half-filled Haldane Chern insulator (chi = 400, Ly = 3)

=========   ========   ========   ========
machine     phi_flow   Ly_flow    V_flow
=========   ========   ========   ========
laptop      --:--:--   --:--:--   --:--:--
---------   --------   --------   --------
Hydra       00:00:00   00:00:00   00:00:00
---------   --------   --------   --------
Piz Daint   00:00:00   00:00:00   00:00:00
=========   ========   ========   ========

3) third-filled Haldane fractional Chern insulator (chi = 500, Ly = 6)

=========   ========   ========   ========
machine     phi_flow   Ly_flow    V_flow
=========   ========   ========   ========
laptop      --:--:--   --:--:--   --:--:--
---------   --------   --------   --------
Hydra       --:--:--   --:--:--   --:--:--
---------   --------   --------   --------
Piz Daint   00:00:00   00:00:00   00:00:00
=========   ========   ========   ========

----

The parameters given for the flows are as follows:

**phi_flow** : phi_min=0, phi_max=1, phi_samp=41

**Ly_flow** : Ly_min=3, Ly_max=6, Ly_samp=2

**V_flow** : V_min=0, V_max=4, V_samp=27