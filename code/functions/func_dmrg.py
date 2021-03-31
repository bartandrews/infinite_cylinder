# --- python imports
import os
import pickle
import gzip
import time
import h5py
import string
# --- TeNPy imports
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
# from tenpy.algorithms.mps_sweeps import OneSiteH, TwoSiteH
from tenpy.tools import hdf5_io
# --- infinite_cylinder imports
import functions.func_proc as fp
from models.heisenberg.heisenberg import HeisenbergModel
from models.ssh.ssh import SSHModel
from models.bbh.bbh import BBHModel
from models.bbh.bbh2 import BBH2Model
from models.haldane.squ_C1 import HalSquC1Model
from models.haldane.hex_C1 import HalHexC1Model
from models.haldane.squ_C2 import HalSquC2Model
from models.haldane.tri_C3 import HalTriC3Model
from models.haldane.squ_CN import HalSquCNModel
from models.hofstadter.squ_1 import HofSqu1Model
from models.hofstadter.hex_1 import HofHex1Model
from models.hofstadter.hex_1_hex_5 import HofHex1Hex5Model
from models.hofstadter.hex_1_hex_5_orbital import HofHex1Hex5OrbitalModel
from models.old.magnetic_lattice.hex_1_hex_5_orbital import FermionicHex1Hex5OrbitalModel  # old code


####################################################
# __get_custom_state (define a custom initial psi) #
####################################################

def __get_custom_state():

    # state = ['up', 'down']*16  # Heisenberg

    # N/4 must be even! e.g. N=8,16,24,32,40,48,etc.

    # state = [1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]*3  # BBH
    # state = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0] * 3  # BBH2
    # state = [0, 1, 0, 1] * int(32 / 4)  # subregion starts at site 1 (SSH)
    state = [0, 1, 1, 0] * int(32 / 4)  # subregion starts at site 2

    return state


####################################################
# __get_product_state (calculates the initial psi) #
####################################################


def __get_product_state(model, ham_params, filling_scale_factor=1, orbital_preference=None):

    nn, nd, LxMUC, Ly, C = ham_params['n'][0], ham_params['n'][1], ham_params['LxMUC'], ham_params['Ly'], ham_params['C']

    # check for Hofstadter model
    if "Hof" in model:
        q = ham_params['nphi'][1]
    else:
        q = 1

    # calculate number of lattice sites
    system_size = int(q) * int(LxMUC) * int(Ly) * int(filling_scale_factor)
    if "Hex" in model or "HalSquC1" in model:  # 2 lattice sites per unit cell
        system_size = system_size * 2

    # n corresponds to the filling fraction of lattice sites in the system
    numb_particles = system_size * int(nn) / int(nd)

    # enforce integer number of particles
    if numb_particles.is_integer():
        numb_particles = int(numb_particles)
    else:
        raise ValueError("Cannot fit an integer number of particles into lattice geometry.")

    # number of sites taken up by a particle (floored fraction)
    numb_sites_per_particle = int(system_size / numb_particles)

    # define empty and full
    empty = '0' if "Bos" in model else 'empty'
    full = '1' if "Bos" in model else 'full'

    if "Orbital" in model or model.endswith("HalSquC2") or model.endswith("HalTriC3"):  # 2 orbitals per site
        # define orbitals
        if "Orbital" in model:
            o1, o2 = 'x', 'y'
        elif model.endswith("HalSquC2") or model.endswith("HalTriC3"):
            o1, o2 = 'A', 'B'
        # define empty_site
        empty_site = [f'{empty}_{o1} {empty}_{o2}']
        if orbital_preference in [f'polarized{o1}', None]:  # default
            lattice_site = [f'{full}_{o1} {empty}_{o2}']
            state = (lattice_site + empty_site * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == f'polarized{o2}':
            lattice_site = [f'{empty}_{o1} {full}_{o2}']
            state = (lattice_site + empty_site * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == 'unpolarized':
            lattice_site_1 = [f'{full}_{o1} {empty}_{o2}']
            lattice_site_2 = [f'{empty}_{o1} {full}_{o2}']
            state = (lattice_site_1 + empty_site * (numb_sites_per_particle - 1)
                     + lattice_site_2 + empty_site * (numb_sites_per_particle - 1)) \
                     * int((system_size / numb_sites_per_particle))
            state = state[:system_size]
        elif orbital_preference == 'filled':
            lattice_site = [f'{full}_{o1} {full}_{o2}']
            state = (lattice_site + empty_site * (2 * numb_sites_per_particle - 1)) \
                     * int((system_size / numb_sites_per_particle))
            state = state[:system_size]
        else:
            raise ValueError("Unknown orbital_preference parameter.")
    elif model.endswith("HalSquCN"):  # N orbitals per site
        # define orbitals
        o = []
        alphabet = list(string.ascii_uppercase)
        for i in range(C):
            o += alphabet[i]
        # define empty_site
        empty_site = f'{empty}_{o[0]}'
        for i in range(1, len(o)):
            empty_site += f' {empty}_{o[i]}'
        empty_site = [empty_site]
        print("empty_site = ", empty_site)

        if orbital_preference in [f'polarized{o[0]}', None]:  # default
            # define lattice_site
            lattice_site = f'{full}_{o[0]}'
            for i in range(1, len(o)):
                lattice_site += f' {empty}_{o[i]}'
            lattice_site = [lattice_site]
            print("lattice_site = ", lattice_site)
            state = (lattice_site + empty_site * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == f'polarized{o[1]}':
            # define lattice_site
            lattice_site = f'{empty}_{o[0]} {full}_{o[1]}'
            for i in range(2, len(o)):
                lattice_site += f' {empty}_{o[i]}'
            lattice_site = [lattice_site]
            print("lattice_site = ", lattice_site)
            state = (lattice_site + empty_site * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == 'unpolarized':
            lattice_sites = []
            for i in range(len(o)):
                lattice_site = ''
                for j in range(len(o)):
                    if j == 0:  # don't prepend space
                        lattice_site += f'{full}_{o[j]}' if i == j else f'{empty}_{o[j]}'
                    else:  # prepend space
                        lattice_site += f' {full}_{o[j]}' if i == j else f' {empty}_{o[j]}'
                lattice_site = [[lattice_site]]
                lattice_sites += lattice_site
            print("lattice_sites = ", lattice_sites)
            state_internal = []
            for i in range(len(o)):
                state_internal += lattice_sites[i] + empty_site * (numb_sites_per_particle - 1)
            state = state_internal * int((system_size / numb_sites_per_particle))
            state = state[:system_size]
        elif orbital_preference == 'filled':
            if numb_particles % C != 0:
                raise ValueError("numb_particles must be a multiple of C for filled orbital_preference.")
            lattice_site = f'{full}_{o[0]}'
            for i in range(1, len(o)):
                lattice_site += f' {full}_{o[i]}'
            lattice_site = [lattice_site]
            print("lattice_site = ", lattice_site)
            state = (lattice_site + empty_site * (2 * numb_sites_per_particle - 1)) \
                    * int((system_size / numb_sites_per_particle))
            state = state[:system_size]
        else:
            raise ValueError("Unknown orbital_preference parameter.")
    else:  # 1 orbital per site
        empty_site = [0]
        if nn / nd <= 0.5:
            state = ([1] + empty_site * (numb_sites_per_particle - 1)) * int((system_size / numb_sites_per_particle))
        else:
            state = [1] * numb_particles

    state += (system_size-len(state)) * empty_site  # top up remainder with empty sites, if necessary

    print("initial state = ", state)
    print("number of particles = ", numb_particles)
    print("number of lattice sites = ", len(state))
    if "Orbital" in model or model.endswith("HalSquC2") or model.endswith("HalTriC3"):  # 2 orbitals per site
        print("number of orbital sites = ", 2 * len(state))
    elif model.endswith("HalSquCN"):  # N orbitals per site
        print("number of orbital sites = ", C * len(state))

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

    if model.endswith("Heisenberg"):
        model_params.clear()
        model_params.update(J=ham_params['t1'])
        model_params.update(D=ham_params['t2'])
        model_params.update(L=ham_params['LxMUC'])
        M = HeisenbergModel(model_params)
    elif model.endswith("SSH"):
        model_params.clear()
        model_params.update(t1=ham_params['t1'])
        model_params.update(t2=ham_params['t2'])
        model_params.update(n=ham_params['n'])
        model_params.update(L=ham_params['LxMUC'])
        M = SSHModel(model_params)
    elif model.endswith("BBH"):
        del model_params['phi'], model_params['mu'], model_params['nphi']
        model_params.update(t2=ham_params['t2'])
        model_params.update(t3=ham_params['t3'])
        model_params.update(t4=ham_params['t4'])
        M = BBHModel(model_params)
    elif model.endswith("BBH2"):
        del model_params['phi'], model_params['mu'], model_params['nphi']
        model_params.update(t2=ham_params['t2'])
        model_params.update(t3=ham_params['t3'])
        model_params.update(t4=ham_params['t4'])
        M = BBH2Model(model_params)
    elif model.endswith("HalSquC1"):
        del model_params['nphi']
        M = HalSquC1Model(model_params)
    elif model.endswith("HalHexC1"):
        del model_params['nphi']
        M = HalHexC1Model(model_params)
    elif model.endswith("HalSquC2"):
        del model_params['nphi']
        M = HalSquC2Model(model_params)
    elif model.endswith("HalTriC3"):
        del model_params['nphi']
        M = HalTriC3Model(model_params)
    elif model.endswith("HalSquCN"):
        del model_params['nphi']
        model_params.update(C=ham_params['C'])
        M = HalSquCNModel(model_params)
    elif model.endswith("HofSqu1"):
        model_params.update(r=ham_params['r'])
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
        # 'update_env': 1000,
        'max_sweeps': 1000,  # default=1000
        'verbose': 5,  # 5
        'N_sweeps_check': 10,  # default=10
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
    __get_product_state(model="FerHofSqu1", ham_params=dict(C=None, n=(4, 7), nphi=(4, 7), LxMUC=1, Ly=3),
                         filling_scale_factor=1, orbital_preference=None)
