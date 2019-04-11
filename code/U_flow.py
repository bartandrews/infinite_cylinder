import numpy as np
import time
import sys
import importlib

import functions as f

parameters_module = "parameters.param_" + str(sys.argv[1])
p = importlib.import_module(parameters_module)


def my_U_flow(model, lattice, initial_state, tile_unit, chi_max, t, mu, V, Lx, Ly, U_min, U_max, U_samp):

    corr_len_stem = f.file_name_stem("corr_len", model, lattice, initial_state, tile_unit, chi_max)
    double_occ_stem = f.file_name_stem("double_occ", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_%s_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat" % (t, U_min, U_max, U_samp, mu, V, Lx, Ly))
    corr_len_file = "data/corr_len/" + corr_len_stem.replace(" ", "_") + leaf
    double_occ_file = "data/double_occ/" + double_occ_stem.replace(" ", "_") + leaf
    open(corr_len_file, "w")
    open(double_occ_file, "w")
    corr_len_data = open(corr_len_file, "a", buffering=1)
    double_occ_data = open(double_occ_file, "a", buffering=1)

    ##################################################################################################################

    engine = f.define_iDMRG_engine_pickle("U_flow", model, lattice, initial_state, tile_unit, chi_max,
                                          t, U_min, mu, V, Lx, Ly, p.use_pickle, p.make_pickle)

    for U in np.linspace(U_min, U_max, U_samp):

        if U != U_min:
            M = f.define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly)
            engine.init_env(model=M)
        engine.run()

        ############
        # corr_len #
        ############

        xi = engine.psi.correlation_length()

        print("{U:.15f}\t{xi:.15f}".format(U=U, xi=xi))
        corr_len_data.write("%.15f\t%.15f\n" % (U, xi))

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
