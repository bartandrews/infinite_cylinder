import time

import functions as f
import parameters as p


def my_ent_spec_real(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors):

    stem = f.file_name_stem("ent_spec_real", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    (E, psi, M) = f.run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

    spectrum = psi.entanglement_spectrum(by_charge=charge_sectors)

    if charge_sectors:

        print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=len(spectrum[0][0][0])))

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        for bond in range(0, Lx*Ly):
            for sector in range(0, len(spectrum[bond])):
                for i in range(0, len(spectrum[bond][sector][1])):
                    print("{charge:d}\t{bond:d}\t{spectrum:.15f}".format(charge=spectrum[bond][sector][0][0], bond=bond,
                                                                         spectrum=spectrum[bond][sector][1][i]))
                    data.write("%i\t%i\t%.15f\n" % (spectrum[bond][sector][0][0],
                                                          bond, spectrum[bond][sector][1][i]))

    else:

        for bond in range(0, Lx * Ly):
            for i in range(0, len(spectrum[bond])):
                print("{bond:d}\t{spectrum:.15f}".format(bond=bond, spectrum=spectrum[bond][i]))
                data.write("%i\t%.15f\n" % (bond, spectrum[bond][i]))


if __name__ == '__main__':

    t0 = time.time()

    my_ent_spec_real(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
                     charge_sectors=True)

    print(time.time() - t0)
