import numpy as np
import time

import functions as f
import parameters as p


def my_ent_spec_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp,
                     charge_sectors):

    stem = f.file_name_stem("ent_spec_flow", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s_%s_%s.dat" % (t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    if charge_sectors:

        engine = f.define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min)

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

            if phi_ext != phi_min:
                M = f.define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)
                engine.init_env(model=M)
            engine.run()

            spectrum = engine.psi.entanglement_spectrum(by_charge=charge_sectors)

            for sector in range(0, len(spectrum[0])):
                for i in range(0, len(spectrum[0][sector][1])):
                    print("{charge:d}\t{phi_ext:.15f}\t{spectrum:.15f}".format(charge=spectrum[0][sector][0][0],
                                                                               phi_ext=phi_ext,
                                                                               spectrum=spectrum[0][sector][1][i]))
                    data.write("%i\t%.15f\t%.15f\n" % (spectrum[0][sector][0][0], phi_ext, spectrum[0][sector][1][i]))

    else:

        engine = f.define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min)

        for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

            M = f.define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, 2*np.pi*phi_ext)
            engine.init_env(model=M)
            engine.run()

            spectrum = engine.psi.entanglement_spectrum(by_charge=charge_sectors)

            for i in range(0, len(spectrum[0])):
                print("{phi_ext:.15f}\t{spectrum:.15f}".format(phi_ext=phi_ext, spectrum=spectrum[0][i]))
                data.write("%.15f\t%.15f\n" % (phi_ext, spectrum[0][i]))


if __name__ == '__main__':

    t0 = time.time()

    my_ent_spec_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
                     phi_min=0, phi_max=1, phi_samp=21, charge_sectors=True)

    print(time.time() - t0)
