import numpy as np
import time

import functions as f
import parameters as p


def my_charge_pump(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp):

    stem = f.file_name_stem("charge_pump", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s_%s_%s.dat" % (t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    engine = f.define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min)

    for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

        if phi_ext != phi_min:
            M = f.define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)
            engine.init_env(model=M)
        engine.run()

        QL_bar = engine.psi.average_charge(bond=0)[0]
        QL = QL_bar  # - engine.psi.get_total_charge()[0]

        print("{phi_ext:.15f}\t{QL:.15f}".format(phi_ext=phi_ext, QL=QL))
        data.write("%.15f\t%.15f\n" % (phi_ext, QL))


if __name__ == '__main__':

    t0 = time.time()

    my_charge_pump(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
                   phi_min=0, phi_max=1, phi_samp=41)

    print(time.time() - t0)
