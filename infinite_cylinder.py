import numpy as np
import random
import sys
import time

from tenpy.networks.mps import MPS
from tenpy.models.fermions_hubbard import FermionicHubbardModel
from models.fermions_haldane import FermionicHaldaneModel
from models.bosons_haldane import BosonicHaldaneModel
from models.fermions_TBG1 import FermionicTBG1Model
from models.fermions_TBG2 import FermionicTBG2Model
from tenpy.algorithms import dmrg


def file_name_stem(tool, model, lattice, initial_state, tile_unit, chi_max, charge_sectors=False):

    if model not in ['Hubbard', 'BosonicHaldane', 'FermionicHaldane', 'TBG1', 'TBG2']:
        sys.exit('Error: Unknown model.')

    if charge_sectors:
        charge_label = '_charge'
    else:
        charge_label = ''

    stem = ("data/%s/%s%s_%s_%s_%s_tile_%s_%s_chi_%s_"
            % (tool, tool, charge_label, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max))

    return stem


def select_initial_psi(model, lattice, initial_state, tile_unit):

    if lattice == "Square":
        lat_basis = 1
    elif lattice == "Honeycomb":
        lat_basis = 2
    else:
        sys.exit('Error: Unknown lattice.')

    Lx, Ly = model.lat.Ls[0], model.lat.Ls[1]

    if initial_state == 0:
        product_state = [tile_unit[0]] * model.lat.N_sites
    elif initial_state == 1:
        product_state = [tile_unit[1]] * model.lat.N_sites
    elif initial_state == 'random':
        product_state = []
        for i in range(model.lat.N_sites):
            product_state.append(random.choice([tile_unit[0], tile_unit[1]]))
    elif initial_state == 'neel':
        product_state = [tile_unit[0]]
        for i in range(0, lat_basis * Lx * Ly - 1):
            if i % 2 == 0:
                product_state.append(tile_unit[1])
            else:
                product_state.append(tile_unit[0])
    elif initial_state == 'third':
        product_state = [tile_unit[0]]
        for i in range(0, lat_basis * Lx * Ly - 1):
            if i % 3 == 0:
                product_state.append(tile_unit[1])
            else:
                product_state.append(tile_unit[0])
    elif initial_state == 'custom':
        product_state = ['down', 'up', 'down', 'up']
    else:
        sys.exit('Error: Unknown initial_state.')

    return product_state


def define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext=0):

    if model == 'Hubbard':
        model_params = dict(cons_N='N', cons_Sz='Sz', t=t, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicHubbardModel(model_params)

    elif model == 'BosonicHaldane':
        model_params = dict(conserve='N', filling=1/2, t=t, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0, phi_ext=phi_ext)
        M = BosonicHaldaneModel(model_params)

    elif model == 'FermionicHaldane':
        model_params = dict(conserve='N', filling=1/2, t=t, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0, phi_ext=phi_ext)
        M = FermionicHaldaneModel(model_params)

    elif model == 'TBG1':
        model_params = dict(cons_N='N', cons_Sz='Sz', t=t, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG1Model(model_params)

    elif model == 'TBG2':
        model_params = dict(conserve='N', filling=1/2, t=t, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG2Model(model_params)

    return M


def define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_ext=0):

    M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)

    product_state = select_initial_psi(M, lattice, initial_state, tile_unit)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,
        'chi_list': {0: 19, 20: 100, 40: chi_max},
        'trunc_params': {
            'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True,
        #     'N_cache': 40
        # },
        'max_E_err': 1.e-10,
        'verbose': 1,
        'N_sweeps_check': 10
    }

    engine = dmrg.EngineCombine(psi, M, dmrg_params)

    return engine


def run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_ext=0):

    M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)

    product_state = select_initial_psi(M, lattice, initial_state, tile_unit)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,
        'chi_list': {0: 19, 20: 100, 40: chi_max},
        'trunc_params': {
            'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True,
        #     'N_cache': 40
        # },
        'max_E_err': 1.e-10,
        'verbose': 0,
        'N_sweeps_check': 10
    }

    info = dmrg.run(psi, M, dmrg_params)

    E = info['E']

    return E, psi, M


def my_corr_len(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min, V_max, V_samp):

    stem = file_name_stem("corr_len", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_%s_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V_min, V_max, V_samp, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    for V in np.linspace(V_min, V_max, V_samp):

        (E, psi, M) = run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

        xi = psi.correlation_length()

        print("{V:.15f}\t{xi:.15f}".format(V=V, xi=xi))
        data.write("%.15f\t%.15f\n" % (V, xi))


def my_charge_pump(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp):

    stem = file_name_stem("charge_pump", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s_%s_%s.dat" % (t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    engine = define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min)

    for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

        if phi_ext != phi_min:
            M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)
            engine.init_env(model=M)
        engine.run()

        QL_bar = engine.psi.average_charge(bond=0)[0]
        QL = QL_bar  # - engine.psi.get_total_charge()[0]

        print("{phi_ext:.15f}\t{QL:.15f}".format(phi_ext=phi_ext, QL=QL))
        data.write("%.15f\t%.15f\n" % (phi_ext, QL))


def my_ent_scal(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly_min, Ly_max, Ly_samp):

    stem = file_name_stem("ent_scal", model, lattice, initial_state, tile_unit, chi_max)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_%s.dat" % (t, U, mu, V, Lx, Ly_min, Ly_max))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    for Ly in np.linspace(Ly_min, Ly_max, Ly_samp):

        (E, psi, M) = run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

        print("{Ly:.15f}\t{SvN:.15f}\t{Sinf:.15f}".format(Ly=Ly, SvN=psi.entanglement_entropy()[0],
                                                       Sinf=psi.entanglement_entropy(n=np.inf)[0]))
        data.write("%.15f\t%.15f\t%.15f\n" % (Ly, psi.entanglement_entropy()[0],
                                           psi.entanglement_entropy(n=np.inf)[0]))


def my_ent_spec_real(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors):

    stem = file_name_stem("ent_spec_real", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    (E, psi, M) = run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

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

    stem = file_name_stem("ent_spec_mom", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    (E, psi, M) = run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

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

    stem = file_name_stem("ent_spec_flow", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s_%s_%s.dat" % (t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    if charge_sectors:

        engine = define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min)

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

            if phi_ext != phi_min:
                M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)
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

        engine = define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min)

        for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

            M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, 2*np.pi*phi_ext)
            engine.init_env(model=M)
            engine.run()

            spectrum = engine.psi.entanglement_spectrum(by_charge=charge_sectors)

            for i in range(0, len(spectrum[0])):
                print("{phi_ext:.15f}\t{spectrum:.15f}".format(phi_ext=phi_ext, spectrum=spectrum[0][i]))
                data.write("%.15f\t%.15f\n" % (phi_ext, spectrum[0][i]))


def my_ent_spec_V_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min, V_max, V_samp,
                       charge_sectors):

    stem = file_name_stem("ent_spec_V_flow", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_%s_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V_min, V_max, V_samp, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    if charge_sectors:

        engine = define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V_min, Lx, Ly)

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        for V in np.linspace(V_min, V_max, V_samp):

            if V != V_min:
                M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly)
                engine.init_env(model=M)
            engine.run()

            spectrum = engine.psi.entanglement_spectrum(by_charge=charge_sectors)

            for sector in range(0, len(spectrum[0])):
                for i in range(0, len(spectrum[0][sector][1])):
                    print("{charge:d}\t{V:.15f}\t{spectrum:.15f}".format(charge=spectrum[0][sector][0][0],
                                                                         V=V, spectrum=spectrum[0][sector][1][i]))
                    data.write("%i\t%.15f\t%.15f\n" % (spectrum[0][sector][0][0], V, spectrum[0][sector][1][i]))

    else:

        engine = define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V_min, Lx, Ly)

        for V in np.linspace(V_min, V_max, V_samp):

            if V != V_min:
                M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly)
                engine.init_env(model=M)
            engine.run()

            spectrum = engine.psi.entanglement_spectrum(by_charge=charge_sectors)

            for i in range(0, len(spectrum[0])):
                print("{V:.15f}\t{spectrum:.15f}".format(V=V, spectrum=spectrum[0][i]))
                data.write("%.15f\t%.15f\n" % (V, spectrum[0][i]))


if __name__ == '__main__':

    # configuration parameters
    model = 'FermionicHaldane'
    lattice = 'Honeycomb'
    initial_state = 'neel'

    if model == 'BosonicHaldane':
        tile_unit = ['0', '1']
    elif model == 'FermionicHaldane':
        tile_unit = ['empty', 'full']
    elif model == 'Hubbard':
        tile_unit = ['down', 'up']
    elif model == 'TBG1':
        tile_unit = ['down_px up_py', 'down_px up_py']
    elif model == 'TBG2':
        tile_unit = ['empty_px full_py', 'empty_px full_py']

    chi_max = 100
    # Hamiltonian parameters (U=0 for FermionicHaldane)
    t, mu, V = -1, 0, 1

    if model in ['BosonicHaldane', 'FermionicHaldane']:
        U = 0
    elif model in ['Hubbard', 'TBG1', 'TBG2']:
        U = 0

    # unit cell
    Lx, Ly = 1, 3

    ####################################################################################################################

    t0 = time.time()

    my_charge_pump(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min=0, phi_max=1,
                   phi_samp=10)
    # my_ent_spec_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min=0, phi_max=1,
    #                  phi_samp=7, charge_sectors=True)
    # my_ent_spec_mom(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors=True)
    # my_ent_scal(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly_min=3, Ly_max=6, Ly_samp=2)
    # my_corr_len(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min=0.3, V_max=1, V_samp=26)
    # my_ent_spec_V_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min=0, V_max=2,
    #                    V_samp=20, charge_sectors=True)
    #
    # my_ent_spec_real(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors=True)

    print(time.time() - t0)
