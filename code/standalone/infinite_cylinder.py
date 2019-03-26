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


def my_ent_spec_mom(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors):

    stem = f.file_name_stem("ent_spec_mom", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    (E, psi, M) = f.run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

    (Un, W, q, ov, trunc_err) = psi.compute_K(perm=M.lat, canonicalize=1.e-6, verbose=0)

    if np.abs(np.abs(ov)-1) > 0.1:
        print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))
        print('Warning: State is not invariant under the permutation.')
    else:
        print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))

    if charge_sectors:

        print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=q.charges.shape[1]))

        # q.to_qflat()[i][0] --> q.to_qflat()[i][n] for different charge entries
        for i in range(len(W)):
            print("{q:d}\t{K:.15f}\t{epsilon:.15f}".format(q=q.to_qflat()[i][0], K=np.angle(W[i])/np.pi,
                                                           epsilon=-np.log(np.abs(W[i]))))
            data.write("%i\t%.15f\t%.15f\n" % (q.to_qflat()[i][0], np.angle(W[i])/np.pi, -np.log(np.abs(W[i]))))

    else:

        for i in range(len(W)):
            print("{K:.15f}\t{epsilon:.15f}".format(K=np.angle(W[i]) / np.pi, epsilon=-np.log(np.abs(W[i]))))
            data.write("%.15f\t%.15f\n" % (np.angle(W[i]) / np.pi, -np.log(np.abs(W[i]))))


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

    # my_charge_pump(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
    #                phi_min=0, phi_max=1, phi_samp=41)
    # my_ent_spec_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
    #                  phi_min=0, phi_max=1, phi_samp=21, charge_sectors=True)
    my_ent_spec_mom(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
                    charge_sectors=True)
    # my_ent_scal(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx,
    #             Ly_min=3, Ly_max=6, Ly_samp=2)
    # my_ent_spec_V_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.Lx, p.Ly,
    #                    V_min=0, V_max=2, V_samp=21, charge_sectors=True)
    # my_corr_len(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.Lx, p.Ly,
    #             V_min=0.3, V_max=1, V_samp=27)
    # my_ent_spec_real(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
    #                  charge_sectors=True)

    print(time.time() - t0)
