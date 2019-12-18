import numpy as np
import time
import sys
import importlib
import tenpy.tools.process as prc

import functions as f

parameters_module = "parameters.param_" + str(sys.argv[1])
p = importlib.import_module(parameters_module)


def my_invt2dash_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, tag,
                      invt2dash_min, invt2dash_max, invt2dash_samp):

    corr_len_invt2dash_flow_stem = \
        f.file_name_stem("corr_len_invt2dash_flow", model, lattice, initial_state, tile_unit, chi_max)
    ent_spec_invt2dash_flow_stem = \
        f.file_name_stem("ent_spec_invt2dash_flow", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_invt2dash_%s_%s_%s_Lx_%s_Ly_%s.dat.%s" %
            (t, U, mu, V, invt2dash_min, invt2dash_max, invt2dash_samp, Lx, Ly, tag))
    corr_len_invt2dash_flow_file = \
        "data/corr_len_invt2dash_flow/" + model + "/" + corr_len_invt2dash_flow_stem.replace(" ", "_") + leaf
    ent_spec_invt2dash_flow_file = \
        "data/ent_spec_invt2dash_flow/" + model + "/" + ent_spec_invt2dash_flow_stem.replace(" ", "_") + leaf
    open(corr_len_invt2dash_flow_file, "w")
    open(ent_spec_invt2dash_flow_file, "w")
    corr_len_invt2dash_flow_data = open(corr_len_invt2dash_flow_file, "a", buffering=1)
    ent_spec_invt2dash_flow_data = open(ent_spec_invt2dash_flow_file, "a", buffering=1)

    ##################################################################################################################

    engine = f.define_iDMRG_engine_pickle("V_flow", model, lattice, initial_state, tile_unit, chi_max,
                                          t, invt2dash, U, mu, V, Lx, Ly, p.use_pickle, p.make_pickle)

    for invt2dash in np.linspace(invt2dash_min, invt2dash_max, invt2dash_samp):

        if invt2dash != invt2dash_min:
            del engine.DMRG_params['chi_list']
            M = f.define_iDMRG_model(model, lattice, t, invt2dash, U, mu, V, Lx, Ly)
            engine.init_env(model=M)
        engine.run()

        ###########################
        # corr_len_invt2dash_flow #
        ###########################

        xi = engine.psi.correlation_length()

        print("{invt2dash:.15f}\t{xi:.15f}".format(invt2dash=invt2dash, xi=xi))
        corr_len_invt2dash_flow_data.write("%.15f\t%.15f\n" % (V, xi))

        ###########################
        # ent_spec_invt2dash_flow #
        ###########################

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        spectrum = engine.psi.entanglement_spectrum(by_charge=True)

        bond = 0

        for sector in range(0, len(spectrum[bond])):
            for i in range(0, len(spectrum[bond][sector][1])):
                print("{charge:d}\t{invt2dash:.15f}\t{spectrum:.15f}".format(charge=spectrum[bond][sector][0][0],
                                                                             invt2dash=invt2dash,
                                                                             spectrum=spectrum[bond][sector][1][i]))
                ent_spec_invt2dash_flow_data.write("%i\t%.15f\t%.15f\n" % (spectrum[bond][sector][0][0], invt2dash,
                                                                           spectrum[bond][sector][1][i]))


if __name__ == '__main__':

    prc.mkl_set_nthreads(1)

    t0 = time.time()

    my_invt2dash_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
                      p.tag, invt2dash_min=10, invt2dash_max=30, invt2dash_samp=21)

    print(time.time() - t0)
