import numpy as np
import time

import functions as f
import parameters as p


def my_ent_scal(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly_min, Ly_max, Ly_samp):

    stem = f.file_name_stem("ent_scal", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_%s.dat" % (t, U, mu, V, Lx, Ly_min, Ly_max))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    for Ly in np.linspace(Ly_min, Ly_max, Ly_samp):

        (E, psi, M) = f.run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

        print("{Ly:.15f}\t{SvN:.15f}\t{Sinf:.15f}".format(Ly=Ly, SvN=psi.entanglement_entropy()[0],
                                                       Sinf=psi.entanglement_entropy(n=np.inf)[0]))
        data.write("%.15f\t%.15f\t%.15f\n" % (Ly, psi.entanglement_entropy()[0],
                                           psi.entanglement_entropy(n=np.inf)[0]))


if __name__ == '__main__':

    t0 = time.time()

    my_ent_scal(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx,
                Ly_min=3, Ly_max=6, Ly_samp=2)

    print(time.time() - t0)
