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

* func_int = interaction functions

    Set of functions to aid in computing the offsite interaction term.

* func_obser = observables functions

    Functions to compute the observables for a ground state, as well as for defining the scalar and nonscalar grouping.

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
- n (numerator= ``nnvalue`` underscore denominator= ``ndvalue``, only range over denominator currently implemented)
- nphi (numerator= ``pvalue`` underscore denominator= ``qvalue``, only range over denominator currently implemented)
- LxMUC
- Ly
- phi
- tag (optional)

NB: For a range of parameter values in an output file, we denote this by the order: min value _ max value _ number of samples (e.g. ``V_0_1_4_Coulomb_1``).

*name = stem + leaf*

Example:  ``data/charge_pump/BosHofSqu1/charge_pump_BosHofSqu1_chi_50_t1_1_t2_0_t2dash_0_U_0_mu_0_V_0_Coulomb_0_n_1_8_8_1_nphi_1_4_4_1_LxMUC_1_Ly_4_4_1_phi_0_2_21.dat``

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

Algorithm scaling
-----------------

Upper-limit scaling relations (actually slightly better due to matrix multiplication optimizations in LAPACK):

Run time: ~O(chi^3 D d^3 + chi^2 D^2 d^2)

Memory usage: ~O(chi^2 d N + 2 chi^2 D N)

* chi = MPS bond dimension
* D = maximum MPO bond dimension
* d = single-site Hilbert space dimension
* N = total number of sites (including extra_dof sites) in the MPS unit cell

Getting started: Madhav Mohan
-----------------------------

1. **Fork the github repository.** You should fork this repository into the directory ``~/PycharmProjects/``. Guide to forking is here: https://guides.github.com/activities/forking/ Please do not submit pull requests or try to push changes to the repository for now. Further useful commands for git versioning can be found in ``~/PycharmProjects/infinite_cylinder/notes/git_commands``.

2. **Set up conda environment (optional).** If you would like to use the exact same conda environment as me, you can now go to ``~/PycharmProjects/infinite_cylinder/notes/`` and type:

``(base) user@computer:~/PycharmProjects/infinite_cylinder/notes/$ conda create --name Bart --file Bart-spec-file.txt``

After you press enter, you will see:

``(Bart) user@computer:~/PycharmProjects/infinite_cylinder/notes/$``

As you can notice, the environment has now changed from base to Bart. Further useful commands for conda environments can be found in ``~/PycharmProjects/infinite_cylinder/notes/conda_commands``.

3. **Configure the PyCharm project.** Go to ``~/PycharmProjects/infinite_cylinder/`` and type:

``user@computer:~/PycharmProjects/infinite_cylinder/$ pycharm-professional &``

This should start an infinite_cylinder Pycharm project. Go to ``File>Settings>Project Interpreter`` and make sure that you have an anaconda project interpreter selected (either the base or Bart). Go to ``File>Settings>Project Structure`` and mark the ``code`` folder as a source folder (it should be blue), and ``Add Content Root`` then add the path to TeNPy directory (for me it is ``/home/bart/TeNPy/``).

4. **Compute your first ground state wavefunction.** Open ``code/ground_state.py`` and run it with the default parameters. You should understand what they all mean. This should take a few minutes to run. When this is done, use a terminal to navigate to ``~/PycharmProjects/infinite_cylinder/pickles`` and notice that a directory has now been created called ground_state. Inside this directory is the ground state wavefunction that you have produced. It is not human readable.

5. **Compute the entanglement entropy of the ground state wavefunction.** In PyCharm, open ``code/observables.py`` and run it with the default parameters. These parameters must match exactly the parameters that you used when you ran ``ground_state.py`` because the code is now looking for a wavefunction file with those parameters in the name. This should take a few seconds. You should see the von Neumann entanglement entropy output to the console. Congratulations, this is effectively a data point on your graph. You know what S is, and you can calculate Ly/lB based on the parameters you gave. Does this value converge as you increase ``chi_max``?

Masters project: Madhav Mohan
-----------------------------

1. **Reproduce an equivalent of Fig. 3 from [Schoond19] for the FerHofSqu1 model at 1/3 filling.** For this, you should use workflow 2 and for each system with filling nu=n/nphi=1/3: compute the von Neumann entanglement entropy, S, for various MPS bond dimensions, chi. That is for fermions with nearest-neighbor interactions: V=10, Vtype='Coulomb', Vrange=1. What do you notice when you plot S vs. 1/chi ? You should see a convergence of the entanglement entropy as you increase the MPS bond dimension (e.g. chi=50, 100, ..., 500). In each case, extrapolate this convergence to get an estimate (with errors) for S in the chi->infty limit. This will form one data point (with error bars) on your graph of S against Ly/lB. Repeating this for a variety of systems with different Ly or nphi, you should get a straight line confirming the area law of entanglement. The (absolute value of the) y-intercept of this straight line is the topological entanglement entropy. What value do you get for the topological entanglement entropy? For the 1/3 state, this value should be 0.549. Keep improving the data points on this plot until you get an agreement to 2 decimal places.

2. **Plot the area law graph for the BosHofSqu1 model at 1/2 filling.** Reproduce the area law plot, as above, now for the BosHofSqu1 model at 1/2 filling. That is hardcore bosons with V=0, Vtype='Coulomb', Vrange=0. You should notice that the computations are faster than for fermions. The topological entanglement entropy for this system is 0.347. Keep improving the data points on this plot until you get an agreement to 2 decimal places.

3. **Decide on a routine.** You have now computed the area law plots for both bosons and fermions. What difficulties arose during your calculations? How large did you have to make the MPS bond dimension to get a convincing convergence extrapolation of S? In what increments is it most efficient to increase chi? Which values of nphi yield fractional quantum Hall states? Are some flux densities more robust than others? Using everything that you have learned, you need to decide on a routine that you can repeat for other Hamiltonians. This project is all about comparing topological entanglement entropy values. So in order to make it a fair test, we need to systematically produce area law plots to the same precision.

4. **Tune the interaction range for the FerHofSqu1 model at 1/3 filling.** You already have a plot for fermions with NN interactions from step 1. Now you can plot the area law for fermions with up to 2nd-NN interactions i.e. V=10, Vtype='Coulomb', Vrange=2. You can then make the interaction longer range by increasing Vrange=3, 4, ..., 10. How far can you get before iDMRG becomes prohibitively slow? What do you notice about the values of the topological entanglement entropy as you increase the interaction range? Originally, for NN-interacting fermions in step 1, you found that the value was 0.549. This value is expected to stay the same since this is the well-known Laughlin state, which is undoubtedly Abelian. Does it indeed stay the same?

5. **[ORIGINAL RESULT] Tune the interaction range for the FerHofSqu1 model at 2/5 filling.** Unlike the 1/3 Laughlin state, the statistics of the FQH state at 2/5 filling is disputed. Jain's composite fermion theory predicts that this ground state has Abelian statistics, whereas the conformal field theory approach (i.e. the Gaffnian) yields non-Abelian statistics. In the recent paper by [Yang19], they claim that this discrepancy is due to the fact that Jain's theory makes an implicit assumption of short-range interactions. The topological entanglement entropy can tell us whether the statistics are Abelian or non-Abelian. Compute the area law now for the FerHofSqu1 model at 2/5 filling for a variety of interaction ranges. What do you notice about the topological entanglement entropy? If their theory is correct, you should observe that the initial value of 0.549 increases as we increase the range of the interactions. Does it increase? You can compare your short-range results with the paper by [Estienne15].

6. **[ORIGINAL RESULT] Tune the interaction range for the FerHofSqu1 model at 3/7 filling.** Following the future work section in the paper by [Yang19], it is now interesting to investigate another disputed filling factor: 3/7. Repeat the calculation from step 5, with this filling. In this case, it is even an original result at short-range. The topological entanglement entropy in the Abelian case is 0.973. Do you get an agreement with this? What happens to this value as you increase the interaction range? Does the topological entanglement entropy increase, as predicted by [Yang19]?

7. **[ORIGINAL RESULT] Diversify the results.** Now that we have two original investigations (2/5 and 3/7), it is time to collect more data to reinforce our claims. The results can be improved in several ways:

* Perform the calculations with the Yukawa interaction.
* Perform the calculations for the hexagonal Hofstadter model.
* Find contested bosonic FQH states, and then perform the calculations for bosons.

All of the code needed for this project is already implemented. The challenge lies in building enough experience in performing iDMRG calculations to understand when the results can be trusted as we push the algorithm to its limits.

References
----------

[Grushin15] "Characterization and stability of a fermionic ν=1/3 fractional Chern insulator" by Adolfo G. Grushin, Johannes Motruk, Michael P. Zaletel, Frank Pollmann, PRB **91**, 035136 (2015). https://arxiv.org/abs/1407.6985

[Zhu19] "Spin/orbital density wave and Mott insulator in two-orbital Hubbard model on honeycomb lattice" by Zheng Zhu, D. N. Sheng, and Liang Fu, Phys. Rev. Lett. **123**, 087602 (2019). https://arxiv.org/abs/1812.05661

[Schoond19] "Interaction-driven plateau transition between integer and fractional Chern Insulators" by Leon Schoonderwoerd, Frank Pollmann, Gunnar Möller, arXiv pre-print (2019). https://arxiv.org/abs/1908.00988

[Yang19] "Effective Abelian theory from a non-Abelian topological order in ν=2/5 fractional quantum Hall effect" by Bo Yang, Ying-Hai Wu, Zlatko Papic, Phys. Rev. B **100**, 245303 (2019). https://arxiv.org/abs/1907.12572

[Estienne15] "Correlation Lengths and Topological Entanglement Entropies of Unitary and Non-Unitary Fractional Quantum Hall Wavefunctions" by B. Estienne, N. Regnault, B. A. Bernevig, Phys. Rev. Lett. **114**, 186801 (2015). https://arxiv.org/abs/1406.6262
