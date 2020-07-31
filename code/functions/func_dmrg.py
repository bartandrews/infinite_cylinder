# --- python imports
import os
import pickle
import gzip
import time
import h5py
# --- TeNPy imports
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
# from tenpy.algorithms.mps_sweeps import OneSiteH, TwoSiteH
from tenpy.tools import hdf5_io
# --- infinite_cylinder imports
import functions.func_proc as fp
from models.haldane.C_1 import HalC1Model
from models.haldane.C_3 import HalC3Model
from models.hofstadter.squ_1 import HofSqu1Model
from models.hofstadter.hex_1 import HofHex1Model
from models.hofstadter.hex_1_hex_5 import HofHex1Hex5Model
from models.hofstadter.hex_1_hex_5_orbital import HofHex1Hex5OrbitalModel
from models.old.magnetic_lattice.hex_1_hex_5_orbital import FermionicHex1Hex5OrbitalModel  # old code


####################################################
# __get_custom_state (define a custom initial psi) #
####################################################

def __get_custom_state():

    state = ['1_A 0_B', '1_A 0_B', '1_A 0_B', '1_A 0_B', '1_A 0_B', '1_A 0_B', '0_A 0_B', '0_A 0_B', '0_A 0_B', '0_A 0_B', '0_A 0_B', '0_A 0_B']

    return state


####################################################
# __get_product_state (calculates the initial psi) #
####################################################


def __get_product_state(model, ham_params, filling_scale_factor=1, orbital_preference=None):

    nn, nd, LxMUC, Ly = ham_params['n'][0], ham_params['n'][1], ham_params['LxMUC'], ham_params['Ly']
    if "Hof" in model:
        q = ham_params['nphi'][1]
    elif model.endswith("HalC3"):
        q = 3
    else:
        q = 1

    if "Orbital" in model:
        numb_particles = 2 * int(q) * int(LxMUC) * int(Ly) * int(filling_scale_factor) * int(nn) / int(nd)
    else:
        numb_particles = int(q) * int(LxMUC) * int(Ly) * int(filling_scale_factor) * int(nn) / int(nd)

    if not numb_particles.is_integer():
        raise ValueError("Cannot fit an integer number of particles into lattice geometry.")
    else:
        numb_particles = int(numb_particles)

    if "Squ" in model or model.endswith("HalC3"):
        system_size = int(q * LxMUC * Ly)
    else:  # "Hex" or "HalC1"
        system_size = int(2 * q * LxMUC * Ly)

    numb_sites_per_particle = int(system_size / numb_particles)

    if "Orbital" in model:
        empty_site = ['0_x 0_y'] if "Bos" in model else ['empty_x empty_y']
        if orbital_preference in ['polarizedx', None]:  # default
            lattice_site = ['1_x 0_y'] if "Bos" in model else ['full_x empty_y']
            state = (lattice_site + empty_site * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == 'polarizedy':
            lattice_site = ['0_x 1_y'] if "Bos" in model else ['empty_x full_y']
            state = (lattice_site + empty_site * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == 'unpolarized':
            if "Bos" in model:
                lattice_site_1 = ['1_x 0_y']
                lattice_site_2 = ['0_x 1_y']
            else:
                lattice_site_1 = ['full_x empty_y']
                lattice_site_2 = ['empty_x full_y']
            state = (lattice_site_1 + empty_site * (numb_sites_per_particle - 1)
                     + lattice_site_2 + empty_site * (numb_sites_per_particle - 1)) \
                     * int((system_size / numb_sites_per_particle))
            state = state[:system_size]
        elif orbital_preference == 'filled':
            lattice_site = ['1_x 1_y'] if "Bos" in model else ['full_x full_y']
            state = (lattice_site + empty_site * (2 * numb_sites_per_particle - 1)) \
                     * int((system_size / numb_sites_per_particle))
            state = state[:system_size]
        else:
            raise ValueError("Unknown orbital_preference parameter.")
    elif model.endswith("HalC3"):
        empty_site = ['0_A 0_B'] if "Bos" in model else ['empty_A empty_B']
        if orbital_preference in ['polarizedA', None]:  # default
            lattice_site = ['1_A 0_B'] if "Bos" in model else ['full_A empty_B']
            state = (lattice_site + empty_site * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == 'polarizedB':
            lattice_site = ['0_A 1_B'] if "Bos" in model else ['empty_A full_B']
            state = (lattice_site + empty_site * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == 'unpolarized':  # old default
            if "Bos" in model:
                lattice_site_1 = ['1_A 0_B']
                lattice_site_2 = ['0_A 1_B']
            else:
                lattice_site_1 = ['full_A empty_B']
                lattice_site_2 = ['empty_A full_B']
            state = (lattice_site_1 + empty_site * (numb_sites_per_particle - 1)
                     + lattice_site_2 + empty_site * (numb_sites_per_particle - 1)) \
                    * int((system_size / numb_sites_per_particle))
            state = state[:system_size]
        elif orbital_preference == 'filled':
            lattice_site = ['1_A 1_B'] if "Bos" in model else ['full_A full_B']
            state = (lattice_site + empty_site * (2 * numb_sites_per_particle - 1)) \
                    * int((system_size / numb_sites_per_particle))
            state = state[:system_size]
        else:
            raise ValueError("Unknown orbital_preference parameter.")
    else:
        empty_site = [0]
        state = ([1] + empty_site * (numb_sites_per_particle - 1)) * int((system_size / numb_sites_per_particle))

    state += (system_size-len(state)) * empty_site  # top up remainder with empty sites, if necessary

    print("initial state = ", state)
    print("number of particles = ", numb_particles)
    print("number of lattice sites = ", len(state))
    if "Orbital" in model or model.endswith("HalC3"):
        print("number of orbital sites = ", 2 * len(state))

    return state


###################################################################
# define_iDMRG_model (defines the parameters for the iDMRG model) #
###################################################################


def define_iDMRG_model(model, ham_params):
    model_params = dict(conserve='N', bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle', verbose=1)
    basic_params = {k: ham_params[k] for k in ('Nmax', 't1', 'mu', 'n', 'nphi',  'LxMUC', 'Ly')}
    if ham_params['V'] != 0:
        basic_params.update({k: ham_params[k] for k in ('V', 'Vtype', 'Vrange')})
    if 'phi' in ham_params:
        basic_params.update({'phi': ham_params['phi']})
    model_params.update(basic_params)

    if "Bos" in model:
        model_params.update(statistics='bosons')
    else:  # "Fer"
        model_params.update(statistics='fermions')

    if model.endswith("HalC1"):
        del model_params['nphi']
        M = HalC1Model(model_params)
    elif model.endswith("HalC3"):
        del model_params['nphi']
        M = HalC3Model(model_params)
    elif model.endswith("HofSqu1"):
        M = HofSqu1Model(model_params)
    elif model.endswith("HofHex1"):
        M = HofHex1Model(model_params)
    elif model.endswith("HofHex1Hex5"):
        model_params.update(t5=ham_params['t5'])
        M = HofHex1Hex5Model(model_params)
    else:  # "HofHex1Hex5Orbital"
        if model == "FerHofHex1Hex5OrbitalOld":  # old code
            del model_params['statistics'], model_params['Vtype'], model_params['Vrange'], model_params['mu']
            model_params.update(t5=ham_params['t5'], t5dash=ham_params['t5dash'], U=ham_params['U'], order='default')
            M = FermionicHex1Hex5OrbitalModel(model_params)
        else:
            model_params.update(t5=ham_params['t5'], t5dash=ham_params['t5dash'], U=ham_params['U'])
            M = HofHex1Hex5OrbitalModel(model_params)

    return M


#################################################################################
# my_iDMRG_pickle (defines an iDMRG engine or run iDMRG and determine pickling) #
#################################################################################


def my_iDMRG_pickle(program, path, model, chi_max, ham_params, use_pickle=False, make_pickle=False, run=True):

    if (use_pickle or make_pickle) and not run:
        raise ValueError("Using/making pickles is only implemented for computing the ground state.")

    state_data = None  # initialize the state_data

    if run:
        pickle_stem = fp.file_name_stem("state", model, chi_max)
        pickle_leaf = fp.file_name_leaf("h5", model, ham_params)
        os.makedirs(os.path.join(path, "pickles", f"{program}", f"{model}", ""), exist_ok=True)
        pickle_path = os.path.join(path, "pickles", f"{program}", f"{model}", pickle_stem + pickle_leaf)

        if use_pickle:
            target_pickle_path = fp.largest_chi_pickle(pickle_path, chi_max)
            target_pickle_file = os.path.split(target_pickle_path)[1]
            t1 = time.time()

            if ".h5" in target_pickle_file:
                ext = "h5"
            elif ".pkl" in target_pickle_file:
                ext = "pkl"
            else:
                raise ValueError("Unknown file extension in target_pickle_path.")

            if ext == "h5":
                file1 = h5py.File(target_pickle_path, 'r')
            else:  # "pkl"
                file1 = (gzip.open if fp.is_gz_file(target_pickle_path) else open)(target_pickle_path, 'rb')

            print(f"Reading {ext} file: ", target_pickle_path)
            print("Writing h5 file: ", pickle_path)
            print(f"Time to open {ext} file (seconds) =", time.time() - t1)
            t2 = time.time()

            if ext == "h5":
                state_data = hdf5_io.load_from_hdf5(file1)
            else:  # "pkl"
                state_data = pickle.load(file1)

            print(f"Time to load {ext} file (seconds) =", time.time() - t2)
            if target_pickle_path != pickle_path.replace(".h5", f".{ext}"):  # improving on a complete pickle
                if state_data['info']['shelve']:
                    raise ValueError(f"Cannot improve on a shelved {ext} file. "
                                     f"Complete the {ext} file first by running at the same chi.")
            else:  # improving on a shelved pickle
                if not state_data['info']['shelve']:
                    raise ValueError(f"Complete {ext} file already exists at this value of chi.")

            file1.close()

    my_iDMRG_output = __my_iDMRG(model, chi_max, ham_params, state_data, run=(True if run else False))

    if make_pickle:
        with h5py.File(pickle_path, 'w') as file2:
            hdf5_io.save_to_hdf5(file2, my_iDMRG_output)

    return my_iDMRG_output  # state_data or engine


#####################################
# __my_iDMRG (implements the iDMRG) #
#####################################


def __my_iDMRG(model, chi_max, ham_params, state_data, run=True):

    if state_data is not None:
        M = state_data['M']
        psi = state_data['psi']
    else:
        M = define_iDMRG_model(model, ham_params)
        orb_pref = 'unpolarized'
        product_state = __get_product_state(model, ham_params) if not ham_params['custom'] else __get_custom_state()
        print(product_state)
        psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        'lanczos_params': {
            # 'reortho': True,
            # 'N_cache': 40
            'cutoff': 1.e-13  # fixes theta=0 error
        },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},  # may need to change first value ~39
        'max_E_err': 1.e-6,
        'max_S_err': 1.e-6,
        # 'norm_tol': 1.e-6,
        # 'norm_tol_iter': 1000,
        'max_sweeps': 1000,
        'verbose': 5,  # 5
        'N_sweeps_check': 10,
        'diag_method': 'default',
        'max_hours': 24*140  # 20 weeks
    }

    if model == "FerHofHex1Hex5OrbitalOld":  # old code
        del dmrg_params['lanczos_params'], dmrg_params['max_sweeps'], dmrg_params['diag_method'], \
            dmrg_params['max_hours']
        dmrg_params.update(mixer_params={'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30}, max_E_err=1.e-10,
                           verbose=1)

    if state_data is not None:
        if state_data['info']['shelve']:
            dmrg_params['sweep_0'] = state_data['info']['sweeps']
        else:
            del dmrg_params['mixer'], dmrg_params['mixer_params'], dmrg_params['chi_list']
            dmrg_params['trunc_params']['chi_max'] = chi_max
        dmrg_params['init_env_data'] = state_data['init_env_data']

    # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
    t3 = time.time()
    engine = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
    print("Time to define engine (seconds) = ", time.time() - t3)
    if run:
        E, psi = engine.run()
        init_env_data = engine.env.get_initialization_data()
        # info = {'dmrg_params': dmrg_params, 'shelve': engine.shelve, 'sweeps': engine.sweeps,
        #         'update_stats': engine.update_stats, 'sweep_stats': engine.sweep_stats}
        info = {'dmrg_params': dmrg_params, 'shelve': engine.shelve, 'sweeps': engine.sweeps,
                'sweep_stats': engine.sweep_stats}
        state_data = {'E': E, 'psi': psi, 'M': M, 'init_env_data': init_env_data, 'info': info}
        return state_data
    else:
        return engine


if __name__ == "__main__":

    # __get_product_state(model="FerHofSqu1", ham_params=dict(n=(4, 45), nphi=(4, 15), LxMUC=1, Ly=6),
    #                     filling_scale_factor=1, orbital_preference=None)
    __get_product_state(model="FerHalC3", ham_params=dict(n=(1, 1), LxMUC=1, Ly=3),
                         filling_scale_factor=1, orbital_preference=None)
