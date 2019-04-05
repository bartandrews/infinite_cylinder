import numpy as np
import time
import sys
import importlib

import functions as f
import tenpy.linalg.np_conserved as npc

parameters_module = "parameters.param_" + str(sys.argv[1])
p = importlib.import_module(parameters_module)


def my_U_flow(model, lattice, initial_state, tile_unit, chi_max, t, mu, V, Lx, Ly, U_min, U_max, U_samp):

    double_occ_stem = f.file_name_stem("double_occ", model, lattice, initial_state, tile_unit, chi_max)
    double_occ_leaf = ("t_%s_U_%s_%s_%s_Lx_%s_Ly_%s.dat" % (t, U_min, U_max, U_samp, Lx, Ly))
    double_occ_file = "data/double_occ/" + double_occ_stem.replace(" ", "_") + double_occ_leaf
    open(double_occ_file, "w")
    double_occ_data = open(double_occ_file, "a", buffering=1)

    ##################################################################################################################

    engine = f.define_iDMRG_engine_pickle("U_flow", model, lattice, initial_state, tile_unit, chi_max,
                                          t, U_min, mu, V, Lx, Ly, p.use_pickle, p.make_pickle)

    for U in np.linspace(U_min, U_max, U_samp):

        if U != U_min:
            M = f.define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly)
            engine.init_env(model=M)
        engine.run()

        ##############
        # double_occ #
        ##############

        term1 = engine.psi.expectation_value('Ntot Ntot')
        term2 = -2*engine.psi.expectation_value('Ntot')

        term1_2 = term1 + term2

        nd_summand = [x + 1 for x in term1_2]
        nd = np.average(nd_summand)

        print("{U:.15f}\t{nd:.15f}".format(U=U, nd=nd))
        double_occ_data.write("%.15f\t%.15f\n" % (U, nd))


if __name__ == '__main__':

    t0 = time.time()

    my_U_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.mu, p.V, p.Lx, p.Ly,
              U_min=1, U_max=10, U_samp=10)

    print(time.time() - t0)
