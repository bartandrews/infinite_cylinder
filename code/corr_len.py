import numpy as np
import time

import functions as f
import parameters as p


def my_corr_len(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min, V_max, V_samp):

    stem = f.file_name_stem("corr_len", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_%s_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V_min, V_max, V_samp, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    for V in np.linspace(V_min, V_max, V_samp):

        (E, psi, M) = f.run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

        xi = psi.correlation_length()

        print("{V:.15f}\t{xi:.15f}".format(V=V, xi=xi))
        data.write("%.15f\t%.15f\n" % (V, xi))


if __name__ == '__main__':

    t0 = time.time()

    my_corr_len(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.Lx, p.Ly,
                V_min=0.3, V_max=1, V_samp=27)

    print(time.time() - t0)
