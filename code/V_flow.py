import numpy as np
import time
import sys
import importlib

import functions as f

parameters_module = "parameters.param_" + str(sys.argv[1])
p = importlib.import_module(parameters_module)


def my_V_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min, V_max, V_samp):

    corr_len_stem = f.file_name_stem("corr_len", model, lattice, initial_state, tile_unit, chi_max)
    ent_spec_V_flow_stem = f.file_name_stem("ent_spec_V_flow", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_%s_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V_min, V_max, V_samp, Lx, Ly))
    corr_len_file = "data/corr_len/" + corr_len_stem.replace(" ", "_") + leaf
    ent_spec_V_flow_file = "data/ent_spec_V_flow/" + ent_spec_V_flow_stem.replace(" ", "_") + leaf
    open(corr_len_file, "w")
    open(ent_spec_V_flow_file, "w")
    corr_len_data = open(corr_len_file, "a", buffering=1)
    ent_spec_V_flow_data = open(ent_spec_V_flow_file, "a", buffering=1)

    ##################################################################################################################

    engine = f.define_iDMRG_engine_pickle("V_flow", model, lattice, initial_state, tile_unit, chi_max,
                                          t, U, mu, V_min, Lx, Ly, p.use_pickle, p.make_pickle)

    for V in np.linspace(V_min, V_max, V_samp):

        if V != V_min:
            del engine.DMRG_params['chi_list']
            M = f.define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly)
            engine.init_env(model=M)
        engine.run()

        ############
        # corr_len #
        ############

        xi = engine.psi.correlation_length()

        print("{V:.15f}\t{xi:.15f}".format(V=V, xi=xi))
        corr_len_data.write("%.15f\t%.15f\n" % (V, xi))

        ###################
        # ent_spec_V_flow #
        ###################

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        spectrum = engine.psi.entanglement_spectrum(by_charge=True)

        bond = 0

        for sector in range(0, len(spectrum[bond])):
            for i in range(0, len(spectrum[bond][sector][1])):
                print("{charge:d}\t{V:.15f}\t{spectrum:.15f}".format(charge=spectrum[bond][sector][0][0],
                                                                     V=V, spectrum=spectrum[bond][sector][1][i]))
                ent_spec_V_flow_data.write("%i\t%.15f\t%.15f\n" % (spectrum[bond][sector][0][0], V,
                                                                   spectrum[bond][sector][1][i]))


if __name__ == '__main__':

    t0 = time.time()

    my_V_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.Lx, p.Ly,
              V_min=0, V_max=4, V_samp=27)

    print(time.time() - t0)
