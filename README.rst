infinite_cylinder
=================

This code is an experimental set of tools for TeNPy, and in due course these are added to the official TeNPy repository. The code focuses on 2D DMRG models on an infinite cylinder.

Prerequisites: TeNPy 0.5+, gnuplot, python 3.6+

Workflow 1 - Using the `flow` programs
--------------------------------------

Occasionally, when a system is simple or predictable enough, it is possible to run iDMRG and record data at the same time. We call these programs 'flows', since they vary the selected independent variable whilst on-the-fly computing a variety of dependent variables of interest. It is also sometimes advantageous to start with the previous state when performing an adiabatic evolution.

**phi_flow** is a program that smoothly varies the external flux through the cylinder. This is used to identify a topological phase and calculate the Chern number. Since the evolution is adiabatic, this flow reuses the state on each iteration.

**Ly_flow** is a program that varies the circumference length of the cylinder. This is used to calculate topological entropy and CFT edge counting. Since the evolution is not adiabatic, this flow recalculates the state on each iteration.

**U_flow** is a program that varies the on-site interaction strength, as defined in the model Hamiltonian. This is used to characterize phase transitions e.g. metal to insulator. Since the evolution is generally not adiabatic, this flow recalculates the state on each iteration.

**V_flow** is a program that varies the off-site interaction strength, as defined in the model Hamiltonian. This is used to characterize phase transitions e.g. metal to FCI. Since the evolution is generally not adiabatic, this flow recalculates the state on each iteration.

**kappa_flow** is a program that varies a lattice parameter coefficient to tune between Hamiltonians. This is used to characterize phase transitions e.g. Hofstadter to TBG. Since the evolution is generally not adiabatic, this flow recalculates the state on each iteration.

The tools employed for each 'flow' are given in the table below.

==========   =====================
**flow**     **tool**
==========   =====================
phi_flow     * charge_pump
             * ent_spec_flow
             * overlap
----------   ---------------------
Ly_flow      * ent_scal
             * ent_spec_real
             * ent_spec_mom
----------   ---------------------
U_flow       * corr_len_U_flow
             * double_occ_U_flow
----------   ---------------------
V_flow       * corr_len_V_flow
             * ent_spec_V_flow
----------   ---------------------
kappa_flow   * corr_len_kappa_flow
             * ent_spec_kappa_flow
==========   =====================

Workflow 2 - Using the `ground_state` and `observables` programs
----------------------------------------------------------------

In cases where the system is not simple to analyze or unpredictable, we need to instead save the ground state that we produce for each iDMRG run using the `ground_state` program. Afterwards, we can load this state and compute our observables of interest individually, using the `observables` program.

The nonscalar observables that are currently implemented for computation are:

* ent_spec_real
* ent_spec_mom
* corr_func

Functions description
---------------------

* func_dmrg = DMRG functions

    Set of functions to calculate the initial state, define the DMRG model, and execute the DMRG.

* func_obser = observables functions

    Functions to compute the observables for a ground state, as well as for their scalar and nonscalar grouping.

* func_proc = file processing functions

    Set of functions to aid with producing output files.


Tools description
-----------------

* charge_pump = charge pump

    This function is designed to plot the equivalent of Figs. 2.a,c) from [Grushin15].

* ent_spec_flow = entanglement spectrum flow

    This function is designed to plot the equivalent of Figs. 2.b,d) from [Grushin15].

* ent_scal = entanglement scaling

    This function is designed to plot the equivalent of Fig. 3.a) from [Grushin15].

* ent_spec_real = entanglement spectrum in real space

    Entanglement energy as a function of bond in the unit cell.

* ent_spec_mom = entanglement spectrum in momentum space

    This function is designed to plot the equivalent of Fig. 3.b) from [Grushin15].

* corr_len_X_flow = correlation length flow with respect to X

    This function is designed to plot the equivalent of the inset in Fig. 3.c) from [Grushin15].

* ent_spec_X_flow = entanglement spectrum flow with respect to X

    This function is designed to plot the equivalent of Fig. 3.c) from [Grushin15].

* double_occ_X_flow = double occupancy flow with respect to X

    This function is designed to plot the equivalent of Fig. 1 from [Zhu19].

* corr_func = two-particle correlation function

    This function is designed to plot the equivalent of Fig. 6 from [Schoond19].

Models description
------------------

hofstadter/hofstadter contains the parent class for all hofstadter models i.e. lattice models in a perpendicular magnetic field

* hofstadter/squ_1

    Hofstadter model with 1st-NN hoppings on a square lattice

* hofstadter/hex_1

    Hofstadter model with 1st-NN hoppings on a honeycomb lattice

* hofstadter/hex_1_hex_5

    Hofstadter model with 1st- and 5th-NN hoppings on a honeycomb lattice

* hofstadter/hex_1_hex_5_orbital

    Hofstadter model with 1st- and 5th-NN hoppings on a honeycomb lattice and two orbitals per site

Directory structure
-------------------

Below is a description of the directory structure of infinite_cylinder, listed alphabetically.

**.idea** is used to store PyCharm configuration files, in case we would like to make changes to the code using a PyCharm project on a remote computer.

**code** contains the source code, categorized into the several parts. **code/functions** is used to store the auxiliary functions for the main programs. **code/models** is used to store custom MPO Hamiltonian python class files. Basic Hamiltonians are already implemented in TeNPy (e.g. Ising model). However, in this directory we store our own Hamiltonian classes. **code/lattices** is used to store custom lattices python class files. Basic lattices are already implemented in TeNPy (e.g. honeycomb). However, in this directory we store our own lattice classes. **code/standalone** is used to store completely independent scripts that do not require the rest of the TeNPy library to run e.g. band structure calculations, Chern number calculations, and plotting scripts. **code/utilities** is used to store python scripts that are used for debugging or checking models, lattices, or other parts of the main code.

**data** is used to store all of the output dat files, organised into their corresponding tool subdirectories (e.g. **ent_spec_real**). Inside each of the tools subdirectories, you will find the models subdirectories (e.g. **BosHofSqu1**). All necessary directories are created at run-time.

**logs** is used to store all of the stdout and stderr output from each run into their corresponding flows subdirectories (e.g. **phi_flow**). Inside each of the flow subdirectories, you will find the models subdirectories (e.g. **BosHofSqu1**). All necessary directories are created at run-time.

**notes** stores Mathematica notebooks for the analysis of the models, and other miscellaneous text files with memos and ideas for future reference.

**pickles** is used to store all of the saved states and DMRG engines into their corresponding flow subdirectories (e.g. **phi_flow**). Inside each of the flow subdirectories, you will find the models subdirectories (e.g. **BosHofSqu1**). All necessary directories are created at run-time.

**scripts** contains bash and python scripts that are used for processing or plotting output, for example.

File naming convention
----------------------

All output .dat files are named in the following order:

*stem*

- tool (e.g. ``charge_pump``)
- model (e.g. ``BosHofSqu1``)
- chi (i.e. ``chi_max``)
- chi_max_K (only for the ent_spec_mom calculation)

*leaf*

- t1
- t2
- t2dash
- kappa (only for the kappa_flow)
- U
- mu
- V
- Vtype (e.g. ``Coulomb``)
- Vrange (e.g. 2 for interactions up to and including 2nd-NN)
- n (numerator underscore denominator, only range over denominator currently implemented)
- nphi (numerator underscore denominator, only range over denominator currently implemented)
- Lx_MUC
- Ly
- phi
- tag (optional)

NB: For a range of parameter values in an output file, we denote this by the order: min value _ max value _ number of samples (e.g. ``V_0_1_4_Coulomb_1``).

*name = stem + leaf*

Example:  ``data/charge_pump/BosHofSqu1/charge_pump_BosHofSqu1_chi_50_t1_1_t2_0_t2dash_0_U_0_mu_0_V_0_Coulomb_0_n_1_8_8_1_nphi_1_4_4_1_Lx_MUC_1_Ly_4_4_1_phi_0_2_21.dat``

Model naming convention
-----------------------

All models are named as follows:

- particle statistics (e.g. ``Bos``/``Fer`` for bosons/fermions -- 3 letter abbreviation)
- model name (e.g. ``Hof`` for the Hofstadter model -- 3 letter abbreviation)
- hopping terms (e.g. ``Squ1Squ2`` for 1st- and 2nd-NN hoppings on a square lattice -- 3 letter abbreviations for the lattices, listed from short to long range hopping)
- other degrees of freedom (e.g. ``OrbitalSpin`` -- full names, in alphabetical order)

Example: ``model='FerHofHex1Hex5Orbital'``

Furthermore, all models with the same model name are grouped into their own subdirectories in ``code/models``.

NB:  model class names do not have the particle statistics prefix and are additionally suffixed with ``Model``.

Pickling capability
-------------------

The pickling capability is used to save the state, or initial state ``[E, psi, M]`` or ``engine`` for a flow. For example, you can save an (expensive) initial DMRG wavefunction, so that you can perform a variety of calculations with it at a later stage. You can set the boolean parameters ``use_pickle`` (to use a pickled state) or ``make_pickle`` (to pickle a state for later) in the parameter files. By default, all pickling is set to False in the flows.

References
----------

[Grushin15] "Characterization and stability of a fermionic ν=1/3 fractional Chern insulator" by Adolfo G. Grushin, Johannes Motruk, Michael P. Zaletel, Frank Pollmann, PRB **91**, 035136 (2015). https://arxiv.org/abs/1407.6985

[Zhu19] "Spin/orbital density wave and Mott insulator in two-orbital Hubbard model on honeycomb lattice" by Zheng Zhu, D. N. Sheng, and Liang Fu, arXiv pre-print (2019). https://arxiv.org/abs/1812.05661

[Schoond19] "Interaction-driven plateau transition between integer and fractional Chern Insulators" by Leon Schoonderwoerd, Frank Pollmann, Gunnar Möller, arXiv pre-print (2019). https://arxiv.org/abs/1908.00988
