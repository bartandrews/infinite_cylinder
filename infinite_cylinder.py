import numpy as np
from tenpy.networks.mps import MPS
from tenpy.models.fermions_hubbard import FermionicHubbardModel
from models.fermions_haldane import FermionicHaldaneModel
from tenpy.algorithms import dmrg
import random
import sys


def select_initial_psi(model, lattice, initial_state, tile_unit):

    if lattice == "Square":
        lat_basis = 1
    elif lattice == "Honeycomb":
        lat_basis = 2
    else:
        sys.exit('Error: Unknown lattice.')

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


def run_iDMRG_Hubbard(initial_state, lattice, tile_unit, chi_max, t, U, mu, V, Lx, Ly):

    model_params = dict(cons_N='N', cons_Sz='Sz', t=t, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                        order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
    M = FermionicHubbardModel(model_params)

    (E, psi, M) = run_iDMRG(lattice, initial_state, tile_unit, chi_max, M)

    return E, psi, M


def run_iDMRG_Haldane(initial_state, lattice, tile_unit, chi_max, t, mu, V, Lx, Ly, phi_ext=0):

    model_params = dict(conserve='N', filling=1/3, t=t, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                        order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0, phi_ext=phi_ext)
    M = FermionicHaldaneModel(model_params)

    (E, psi, M) = run_iDMRG(lattice, initial_state, tile_unit, chi_max, M)

    return E, psi, M


def run_iDMRG(lattice, initial_state, tile_unit, chi_max, model):

    product_state = select_initial_psi(model, lattice, initial_state, tile_unit)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(model.lat.mps_sites(), product_state, bc=model.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,
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

    info = dmrg.run(psi, model, dmrg_params)

    E = info['E']

    return E, psi, model


def my_corr_len(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min, V_max, V_samp):

    directory = "data/corr_len/"

    if model == 'Hubbard':
        file = ("corr_len_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_U_%s_mu_%s_V_%s_%s_%s_Lx_%s_Ly_%s.dat"
                % (model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, U, mu, V_min, V_max, V_samp,
                   Lx, Ly))
    elif model == 'Haldane':
        file = ("corr_len_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_mu_%s_V_%s_%s_%s_Lx_%s_Ly_%s.dat"
                % (model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, mu, V_min, V_max, V_samp,
                   Lx, Ly))
    else:
        sys.exit('Error: Unknown model.')

    dat_file = directory + file

    open(dat_file, "w")
    data = open(dat_file, "a")

    for V in np.linspace(V_min, V_max, V_samp):

        if model == 'Hubbard':
            (E, psi, M) = run_iDMRG_Hubbard(initial_state, lattice, tile_unit, chi_max, t, U, mu, V, Lx, Ly)
        elif model == 'Haldane':
            (E, psi, M) = run_iDMRG_Haldane(initial_state, lattice, tile_unit, chi_max, t, mu, V, Lx, Ly)
        else:
            sys.exit('Error: Unknown model.')

        numb_sites = len(M.lat.mps_sites())
        xi = psi.correlation_length()

        print("{V:.15f}\t{xi:.15f}".format(V=V, xi=xi))
        data.write("%.15f\t%.15f\n" % (V, xi))


def my_charge_pump(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp):

    directory = "data/charge_pump/"

    if model == 'Hubbard':
        file = ("charge_pump_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s_%s_%s.dat"
                % (model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, U, mu, V,
                   Lx, Ly, phi_min, phi_max, phi_samp))
    elif model == 'Haldane':
        file = ("charge_pump_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s_%s_%s.dat"
                % (model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, mu, V,
                   Lx, Ly, phi_min, phi_max, phi_samp))
    else:
        sys.exit('Error: Unknown model.')

    dat_file = directory + file

    open(dat_file, "w")
    data = open(dat_file, "a")

    for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

        if model == 'Hubbard':
            (E, psi, M) = run_iDMRG_Hubbard(initial_state, lattice, tile_unit, chi_max, t, U, mu, V, Lx, Ly)
        elif model == 'Haldane':
            (E, psi, M) = run_iDMRG_Haldane(initial_state, lattice, tile_unit, chi_max, t, mu, V, Lx, Ly, phi_ext)
        else:
            sys.exit('Error: Unknown model.')

        numb_sites = len(M.lat.mps_sites())
        dN = psi.expectation_value('dN')

        # cumulative charges to the left of bonds
        QL_array = []
        for i in range(len(dN)):
            QL_array.append(sum(dN[0:i]))

        # print(psi.get_SL(numb_sites//2))

        # total charge to the left of middle bond
        QL_bar = QL_array[numb_sites//2]

        # physical total charge to the left of middle bond (for the iMPS)
        QL = QL_bar - np.mean(QL_array)
        print("{phi_ext:.15f}\t{QL:.15f}".format(phi_ext=phi_ext, QL=QL))
        data.write("%.15f\t%.15f\n" % (phi_ext, QL))


def my_ent_scal(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly_min, Ly_max):

    directory = "data/ent_scal/"

    if model == 'Hubbard':
        file = ("ent_scal_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_%s.dat"
                % (model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, U, mu, V,
                   Lx, Ly_min, Ly_max))
    elif model == 'Haldane':
        file = ("ent_scal_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_mu_%s_V_%s_Lx_%s_Ly_%s_%s.dat"
                % (model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, mu, V,
                   Lx, Ly_min, Ly_max))
    else:
        sys.exit('Error: Unknown model.')

    dat_file = directory + file

    open(dat_file, "w")
    data = open(dat_file, "a")

    for Ly_iter in range(Ly_min, Ly_max+1):

        if model == 'Hubbard':
            (E, psi, M) = run_iDMRG_Hubbard(initial_state, lattice, tile_unit, chi_max, t, U, mu, V, Lx, Ly)
        elif model == 'Haldane':
            (E, psi, M) = run_iDMRG_Haldane(initial_state, lattice, tile_unit, chi_max, t, mu, V, lattice, Lx, Ly)
        else:
            sys.exit('Error: Unknown model.')

        print("{Ly:d}\t{SvN:.15f}\t{Sinf:.15f}".format(Ly=Ly_iter, SvN=psi.entanglement_entropy()[(Lx*Ly_iter-1)//2],
                                                       Sinf=psi.entanglement_entropy(n=np.inf)[(Lx*Ly_iter-1)//2]))
        data.write("%i\t%.15f\t%.15f\n" % (Ly_iter, psi.entanglement_entropy()[(Lx*Ly_iter-1)//2],
                                           psi.entanglement_entropy(n=np.inf)[(Lx*Ly_iter-1)//2]))


def my_ent_spec_real(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors):

    if charge_sectors == True:
        charge_label = '_charge'
    elif charge_sectors == False:
        charge_label = ''
    else:
        sys.exit('Error: Unknown charge_sectors.')

    directory = "data/ent_spec_real/"

    if model == 'Hubbard':
        file = ("ent_spec_real%s_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat"
                % (charge_label, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, U, mu, V,
                   Lx, Ly))
        (E, psi, M) = run_iDMRG_Hubbard(initial_state, lattice, tile_unit, chi_max, t, U, mu, V, Lx, Ly)
    elif model == 'Haldane':
        file = ("ent_spec_real%s_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat"
                % (charge_label, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, mu, V,
                   Lx, Ly))
        (E, psi, M) = run_iDMRG_Haldane(initial_state, lattice, tile_unit, chi_max, t, mu, V, Lx, Ly)
    else:
        sys.exit('Error: Unknown model.')

    dat_file = directory + file

    open(dat_file, "w")
    data = open(dat_file, "a")

    spectrum = psi.entanglement_spectrum(by_charge=charge_sectors)

    if charge_sectors == False:

        for bond in range(0, Lx*Ly):
            for i in range(0, len(spectrum[bond])):
                print("{bond:d}\t{spectrum:.15f}".format(bond=bond, spectrum=spectrum[bond][i]))
                data.write("%i\t%.15f\n" % (bond, spectrum[bond][i]))

    elif charge_sectors == True:

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
        sys.exit('Error: Unknown charge_sectors.')


def my_ent_spec_mom(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors):

    if charge_sectors == True:
        charge_label = '_charge'
    elif charge_sectors == False:
        charge_label = ''
    else:
        sys.exit('Error: Unknown charge_sectors.')

    directory = "data/ent_spec_mom/"

    if model == 'Hubbard':
        file = ("ent_spec_mom%s_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat"
                % (charge_label, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, U, mu, V,
                   Lx, Ly))
        (E, psi, M) = run_iDMRG_Hubbard(initial_state, lattice, tile_unit, chi_max, t, U, mu, V, Lx, Ly)
    elif model == 'Haldane':
        file = ("ent_spec_mom%s_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat"
                % (charge_label, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, mu, V,
                   Lx, Ly))
        (E, psi, M) = run_iDMRG_Haldane(initial_state, lattice, tile_unit, chi_max, t, mu, V, Lx, Ly)
    else:
        sys.exit('Error: Unknown model.')

    dat_file = directory + file

    open(dat_file, "w")
    data = open(dat_file, "a")

    (Un, W, q, ov, trunc_err) = psi.compute_K(perm=M.lat, canonicalize=1.e-3, verbose=0)

    if np.abs(np.abs(ov)-1) > 0.1:
        print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))
        print('Warning: State is not invariant under the permutation.')
    else:
        print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))

    if charge_sectors == False:

        for i in range(len(W)):
            print("{K:.15f}\t{epsilon:.15f}".format(K=np.angle(W[i])/np.pi, epsilon=-np.log(np.abs(W[i]))))
            data.write("%.15f\t%.15f\n" % (np.angle(W[i])/np.pi, -np.log(np.abs(W[i]))))

    elif charge_sectors == True:

        print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=q.charges.shape[1]))

        # q.to_qflat()[i][0] --> q.to_qflat()[i][n] for different charge entries
        for i in range(len(W)):
            print("{q:d}\t{K:.15f}\t{epsilon:.15f}".format(q=q.to_qflat()[i][0], K=np.angle(W[i])/np.pi,
                                                           epsilon=-np.log(np.abs(W[i]))))
            data.write("%i\t%.15f\t%.15f\n" % (q.to_qflat()[i][0], np.angle(W[i])/np.pi, -np.log(np.abs(W[i]))))

    else:
        sys.exit('Error: Unknown charge_sectors.')


def my_ent_spec_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min, phi_max, phi_samp,
                     charge_sectors):

    if charge_sectors == True:
        charge_label = '_charge'
    elif charge_sectors == False:
        charge_label = ''
    else:
        sys.exit('Error: Unknown charge_sectors.')

    directory = "data/ent_spec_flow/"

    if model == 'Hubbard':
        file = ("ent_spec_flow%s_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s_%s_%s.dat"
                % (charge_label, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, U, mu, V,
                   Lx, Ly, phi_min, phi_max, phi_samp))
    elif model == 'Haldane':
        file = ("ent_spec_flow%s_%s_%s_%s_tile_%s_%s_chi_%s_t_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s_%s_%s.dat"
                % (charge_label, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max,
                   t, mu, V,
                   Lx, Ly, phi_min, phi_max, phi_samp))
    else:
        sys.exit('Error: Unknown model.')

    dat_file = directory + file

    open(dat_file, "w")
    data = open(dat_file, "a")

    if charge_sectors == False:

        for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

            if model == 'Hubbard':
                (E, psi, M) = run_iDMRG_Hubbard(initial_state, lattice, tile_unit, chi_max, t, U, mu, V, Lx, Ly)
            elif model == 'Haldane':
                (E, psi, M) = run_iDMRG_Haldane(initial_state, lattice, tile_unit, chi_max, t, mu, V, Lx, Ly, phi_ext)
            else:
                sys.exit('Error: Unknown model.')

            spectrum = psi.entanglement_spectrum(by_charge=charge_sectors)

            for i in range(0, len(spectrum[0])):
                print("{phi_ext:.15f}\t{spectrum:.15f}".format(phi_ext=phi_ext, spectrum=spectrum[0][i]))
                data.write("%.15f\t%.15f\n" % (phi_ext, spectrum[0][i]))

    elif charge_sectors == True:

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

            if model == 'Hubbard':
                (E, psi, M) = run_iDMRG_Hubbard(initial_state, lattice, tile_unit, chi_max, t, U, mu, V, Lx, Ly)
            elif model == 'Haldane':
                (E, psi, M) = run_iDMRG_Haldane(initial_state, lattice, tile_unit, chi_max, t, mu, V, Lx, Ly, phi_ext)
            else:
                sys.exit('Error: Unknown model.')

            spectrum = psi.entanglement_spectrum(by_charge=charge_sectors)

            for sector in range(0, len(spectrum[0])):
                for i in range(0, len(spectrum[0][sector][1])):
                    print("{charge:d}\t{phi_ext:.15f}\t{spectrum:.15f}".format(charge=spectrum[0][sector][0][0],
                                                                               phi_ext=phi_ext,
                                                                               spectrum=spectrum[0][sector][1][i]))
                    data.write("%i\t%.15f\t%.15f\n" % (spectrum[0][sector][0][0], phi_ext, spectrum[0][sector][1][i]))

    else:
        sys.exit('Error: Unknown charge_sectors.')


if __name__ == '__main__':

    # configuration parameters
    model = 'Hubbard'
    lattice = 'Square'
    initial_state = 'neel'
    tile_unit = ['up', 'down']
    chi_max = 30
    # Hamiltonian parameters (U only for Hubbard)
    t, U, mu, V = -1, 0, 0, 1
    # unit cell
    Lx, Ly = 2, 2

    # my_corr_len(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, Lx, Ly, V_min=0, V_max=1, V_samp=4)
    my_charge_pump(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min=0, phi_max=1,
                   phi_samp=4)
    # my_ent_scal(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly_min=2, Ly_max=4)
    # my_ent_spec_real(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors=True)
    # my_ent_spec_mom(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors=True)
    # my_ent_spec_flow(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_min=0, phi_max=1,
    #                  phi_samp=4, charge_sectors=True)
