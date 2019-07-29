from tenpy.networks.mps import MPS
from tenpy.models.fermions_hubbard import FermionicHubbardModel
from models.fermions_haldane import FermionicHaldaneModel
from models.fermions_hofstadter import FermionicHofstadterModel
from models.fermions_hex_1 import FermionicHex1Model
from models.fermions_twist import FermionicTwistModel
from models.fermions_complete_twist import FermionicCompleteTwistModel
from models.fermions_pi_flux import FermionicPiFluxModel
from models.fermions_C3_haldane import FermionicC3HaldaneModel
from models.bosons_haldane import BosonicHaldaneModel
from models.bosons_haldane_2 import BosonicHaldane2Model
from models.fermions_TBG1 import FermionicTBG1Model
from models.fermions_TBG2 import FermionicTBG2Model
from models.fermions_TBG3 import FermionicTBG3Model
from models.fermions_TBG4 import FermionicTBG4Model
from models.fermions_TBG5 import FermionicTBG5Model
from models.fermions_TBG6 import FermionicTBG6Model

from tenpy.algorithms import dmrg
# from tenpy.algorithms.mps_sweeps import OneSiteH, TwoSiteH

import sys
import random
import pickle


def file_name_stem(tool, model, lattice, initial_state, tile_unit, chi_max):

    if model not in ['Hubbard', 'BosonicHaldane', 'BosonicHaldane2', 'FermionicHaldane',
                     'FermionicHofstadter', 'FermionicHex1', 'FermionicTwist', 'FermionicCompleteTwist',
                     'FermionicPiFlux', 'FermionicC3Haldane', 'TBG1', 'TBG2', 'TBG3', 'TBG4', 'TBG5', 'TBG6']:
        sys.exit('Error: Unknown model.')

    stem = ("%s_%s_%s_%s_tile_%s_%s_chi_%s_"
            % (tool, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max))

    return stem


def select_initial_psi(model, lattice, initial_state, tile_unit):

    if lattice == "Square" or lattice == "Triangular":
        lat_basis = 1
    elif lattice == "Honeycomb":
        lat_basis = 2
    elif lattice == "TripartiteTriangular":
        lat_basis = 3
    elif lattice == "BipartiteSquare":
        lat_basis = 4
    elif lattice == "FiveBandLattice":
        lat_basis = 6
    elif lattice == "MagneticSquare":
        lat_basis = 3
    elif lattice == "MagneticHoneycomb":
        lat_basis = 10
    elif lattice == "MagneticTwist":
        lat_basis = 14
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
        product_state = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
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

    elif model == 'FermionicHofstadter':
        model_params = dict(conserve='N', t=t, mu=mu, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicHofstadterModel(model_params)

    elif model == 'FermionicHex1':
        model_params = dict(conserve='N', t=t, mu=mu, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicHex1Model(model_params)

    elif model == 'FermionicTwist':
        model_params = dict(conserve='N', t=t, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicTwistModel(model_params)

    elif model == 'FermionicCompleteTwist':
        model_params = dict(conserve='N', t=t, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicCompleteTwistModel(model_params)

    elif model == 'FermionicPiFlux':
        model_params = dict(conserve='N', t=t, mu=mu, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicPiFluxModel(model_params)

    elif model == 'FermionicC3Haldane':
        model_params = dict(conserve='N', t=t, mu=mu, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicC3HaldaneModel(model_params)

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

    elif model == 'TBG4':
        model_params = dict(conserve='N', t=t, U=U, mu=mu, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1)
        M = FermionicTBG4Model(model_params)

    elif model == 'TBG6':
        model_params = dict(conserve='N', t=t, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG6Model(model_params)

    return M


def define_iDMRG_spin_model(model, lattice, J, Js, Jv, Lx, Ly, phi_ext=0):

    if model == 'TBG5':
        model_params = dict(conserve='Sz', J=J, Js=Js, Jv=Jv, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=1, phi_ext=phi_ext)
        M = FermionicTBG5Model(model_params)

    return M


def define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, phi_ext=0):

    M = define_iDMRG_model(model, lattice, t, U, mu, V, Lx, Ly, phi_ext)

    product_state = select_initial_psi(M, lattice, initial_state, tile_unit)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        # 'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True,
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
        'max_E_err': 1.e-10,
        'max_S_err': 1.e-6,
        # 'norm_tol': 1.e-6,
        # 'norm_tol_iter': 1000,
        # 'max_sweeps': 150,
        'verbose': 1,
        'N_sweeps_check': 10
    }

    # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
    engine = dmrg.EngineCombine(psi, M, dmrg_params)

    return engine


def define_iDMRG_spin_engine(model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly, phi_ext=0):

    M = define_iDMRG_spin_model(model, lattice, J, Js, Jv, Lx, Ly, phi_ext)

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
        # 'norm_tol': 1.e-6,
        # 'norm_tol_iter': 1000,
        # 'max_sweeps': 150,
        'verbose': 1,
        'N_sweeps_check': 10
    }

    engine = dmrg.OneSiteDMRGEngine(psi, M, dmrg_params)

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


def define_iDMRG_spin_engine_pickle(flow, model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly,
                                    use_pickle=False, make_pickle=False, phi_ext=0):

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("engine", model, lattice, initial_state, tile_unit, chi_max)
        pickle_leaf = ("J_%s_Js_%s_Jv_%s_Lx_%s_Ly_%s_phi_%s.pkl" % (J, Js, Jv, Lx, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:

        with open(pickle_file, 'rb') as file1:
            engine = pickle.load(file1)

    else:

        engine = define_iDMRG_spin_engine(model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly, phi_ext)

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


def run_iDMRG_spin(model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly, phi_ext=0):

    M = define_iDMRG_spin_model(model, lattice, J, Js, Jv, Lx, Ly, phi_ext)

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


def run_iDMRG_spin_pickle(flow, model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly,
                          use_pickle=False, make_pickle=False, phi_ext=0):

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("E_psi_M", model, lattice, initial_state, tile_unit, chi_max)
        pickle_leaf = ("J_%s_Js_%s_Jv_%s_Lx_%s_Ly_%s_phi_%s.pkl" % (J, Js, Jv, Lx, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:

        with open(pickle_file, 'rb') as file1:
            [E, psi, M] = pickle.load(file1)

    else:

        (E, psi, M) = run_iDMRG_spin(model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly, phi_ext)

        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                pickle.dump([E, psi, M], file2)

    return E, psi, M

