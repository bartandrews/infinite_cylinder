import numpy as np
import time
import sys
import importlib

import functions as f

parameters_module = "parameters.param_" + str(sys.argv[1])
p = importlib.import_module(parameters_module)


def my_J_flow(model, lattice, initial_state, tile_unit, chi_max, J, Jv, Lx, Ly, Js_min, Js_max, Js_samp):

    energy_stem = f.file_name_stem("energy", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("J_%s_J_%s_%s_%s_Jv_%s_Lx_%s_Ly_%s.dat" % (J, Js_min, Js_max, Js_samp, Jv, Lx, Ly))
    energy_file = "data/energy/" + energy_stem.replace(" ", "_") + leaf
    open(energy_file, "w")
    energy_data = open(energy_file, "a", buffering=1)

    ##################################################################################################################

    engine = f.define_iDMRG_spin_engine_pickle("J_flow", model, lattice, initial_state, tile_unit, chi_max,
                                               J, Js_min, Jv, Lx, Ly, p.use_pickle, p.make_pickle)

    for Js in np.linspace(Js_min, Js_max, Js_samp):

        if Js != Js_min:
            del engine.DMRG_params['chi_list']
            M = f.define_iDMRG_spin_model(model, lattice, J, Js, Jv, Lx, Ly)
            engine.init_env(model=M)
        (E, _) = engine.run()

        ##########
        # energy #
        ##########

        energy = E

        print("{Js:.15f}\t{energy:.15f}".format(Js=Js, energy=energy))
        energy_data.write("%.15f\t%.15f\n" % (Js, energy))


if __name__ == '__main__':

    t0 = time.time()

    my_J_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.J, p.Jv, p.Lx, p.Ly,
              Js_min=0.2, Js_max=0.3, Js_samp=2)

    print(time.time() - t0)
