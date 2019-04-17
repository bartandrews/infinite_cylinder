from tenpy.networks.mps import MPS
from tenpy.models.fermions_hubbard import FermionicHubbardModel
from models.fermions_haldane import FermionicHaldaneModel
from models.bosons_haldane import BosonicHaldaneModel
from models.bosons_haldane_2 import BosonicHaldane2Model
from models.fermions_TBG1 import FermionicTBG1Model
from models.fermions_TBG2 import FermionicTBG2Model
from models.fermions_TBG3 import FermionicTBG3Model
from tenpy.algorithms import dmrg

import sys
import random
import pickle


def file_name_stem(tool, model, lattice, initial_state, tile_unit, chi_max):

    if model not in ['Hubbard', 'BosonicHaldane', 'BosonicHaldane2', 'FermionicHaldane', 'TBG1', 'TBG2', 'TBG3']:
        sys.exit('Error: Unknown model.')

    stem = ("%s_%s_%s_%s_tile_%s_%s_chi_%s_"
            % (tool, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max))

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
        model_params = dict(conserve='N', t=t, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0, phi_ext=phi_ext)
        M = BosonicHaldaneModel(model_params)

    elif model == 'BosonicHaldane2':
        model_params = dict(conserve='N', t=t, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0, phi_ext=phi_ext)
        M = BosonicHaldane2Model(model_params)

    elif model == 'FermionicHaldane':
        model_params = dict(conserve='N', t=t, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0, phi_ext=phi_ext)
        M = FermionicHaldaneModel(model_params)

    elif model == 'TBG1':
        model_params = dict(cons_N='N', cons_Sz='Sz', t=t, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG1Model(model_params)

    elif model == 'TBG2':
        model_params = dict(conserve='N', t=t, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG2Model(model_params)

    elif model == 'TBG3':
        model_params = dict(cons_N='N', cons_Sz='Sz', t=t, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG3Model(model_params)

    return M


def define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_ext=0):

    M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)

    product_state = select_initial_psi(M, lattice, initial_state, tile_unit)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True,
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
        'max_E_err': 1.e-10,
        'max_S_err': 1.e-6,
        # 'max_sweeps': 150,
        'verbose': 1,
        'N_sweeps_check': 10
    }

    engine = dmrg.EngineCombine(psi, M, dmrg_params)

    return engine


def define_iDMRG_engine_pickle(flow, model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly,
                               use_pickle=False, make_pickle=False, phi_ext=0):

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("engine", model, lattice, initial_state, tile_unit, chi_max)
        pickle_leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s.pkl" % (t, U, mu, V, Lx, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:

        with open(pickle_file, 'rb') as file1:
            engine = pickle.load(file1)

    else:

        engine = define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_ext)

        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                pickle.dump(engine, file2)

    return engine


def run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_ext=0):

    M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)

    product_state = select_initial_psi(M, lattice, initial_state, tile_unit)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
        'max_E_err': 1.e-10,
        'max_S_err': 1.e-6,
        # 'max_sweeps': 150,
        'verbose': 1,
        'N_sweeps_check': 10
    }

    info = dmrg.run(psi, M, dmrg_params)

    E = info['E']

    return E, psi, M


def run_iDMRG_pickle(flow, model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly,
                     use_pickle=False, make_pickle=False, phi_ext=0):

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("E_psi_M", model, lattice, initial_state, tile_unit, chi_max)
        pickle_leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_phi_%s.pkl" % (t, U, mu, V, Lx, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:

        with open(pickle_file, 'rb') as file1:
            [E, psi, M] = pickle.load(file1)

    else:

        (E, psi, M) = run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_ext)

        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                pickle.dump([E, psi, M], file2)

    return E, psi, M
