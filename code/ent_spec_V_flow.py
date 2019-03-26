import numpy as np
import time

import functions as f
import parameters as p


def my_ent_spec_V_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min, V_max, V_samp,
                       charge_sectors):

    stem = f.file_name_stem("ent_spec_V_flow", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_%s_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V_min, V_max, V_samp, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    if charge_sectors:

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        for V in np.linspace(V_min, V_max, V_samp):

            (E, psi, M) = f.run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

            spectrum = psi.entanglement_spectrum(by_charge=charge_sectors)

            for sector in range(0, len(spectrum[0])):
                for i in range(0, len(spectrum[0][sector][1])):
                    print("{charge:d}\t{V:.15f}\t{spectrum:.15f}".format(charge=spectrum[0][sector][0][0],
                                                                         V=V, spectrum=spectrum[0][sector][1][i]))
                    data.write("%i\t%.15f\t%.15f\n" % (spectrum[0][sector][0][0], V, spectrum[0][sector][1][i]))

    else:

        for V in np.linspace(V_min, V_max, V_samp):

            (E, psi, M) = f.run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

            spectrum = psi.entanglement_spectrum(by_charge=charge_sectors)

            for i in range(0, len(spectrum[0])):
                print("{V:.15f}\t{spectrum:.15f}".format(V=V, spectrum=spectrum[0][i]))
                data.write("%.15f\t%.15f\n" % (V, spectrum[0][i]))


if __name__ == '__main__':

    t0 = time.time()

    my_ent_spec_V_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.Lx, p.Ly,
                       V_min=0, V_max=2, V_samp=21, charge_sectors=True)

    print(time.time() - t0)
