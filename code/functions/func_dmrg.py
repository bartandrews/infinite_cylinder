# --- python imports
import os
import pickle
# --- TeNPy imports
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
# from tenpy.algorithms.mps_sweeps import OneSiteH, TwoSiteH
# --- infinite_cylinder imports
import functions.func_proc as fp
from models.hofstadter.squ_1 import HofSqu1Model
from models.hofstadter.hex_1 import HofHex1Model
from models.hofstadter.hex_1_hex_5 import HofHex1Hex5Model
from models.hofstadter.hex_1_hex_5_orbital import HofHex1Hex5OrbitalModel


####################################################
# __get_product_state (calculates the initial psi) #
####################################################


def __get_product_state(model, ham_params, filling_scale_factor=1, orbital_preference=None):

    nn, nd, q, LxMUC, Ly = \
        ham_params['n'][0], ham_params['n'][1], ham_params['nphi'][1], ham_params['LxMUC'], ham_params['Ly']

    if "Orbital" in model:
        numb_particles = 2 * int(q) * int(LxMUC) * int(Ly) * int(filling_scale_factor) * int(nn) / int(
            nd)
    else:
        numb_particles = int(q) * int(LxMUC) * int(Ly) * int(filling_scale_factor) * int(nn) / int(
            nd)

    if not numb_particles.is_integer():
        raise ValueError("Cannot fit an integer number of particles into lattice geometry.")
    else:
        numb_particles = int(numb_particles)

    if "Squ" in model:
        system_size = int(q * LxMUC * Ly)
    else:  # "Hex"
        system_size = int(2 * q * LxMUC * Ly)

    numb_sites_per_particle = int(system_size / numb_particles)

    if "Orbital" in model:
        empty_site = ['0_x 0_y'] if "Bos" in model else ['empty_x empty_y']
        if orbital_preference in ['polarizedx', None]:
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
                     * int((system_size / (2 * numb_sites_per_particle)))
        elif orbital_preference == 'filled':
            lattice_site = ['1_x 1_y'] if "Bos" in model else ['full_x full_y']
            state = (lattice_site + empty_site * (2 * numb_sites_per_particle - 1)) \
                     * int((system_size / (2 * numb_sites_per_particle)))
        else:
            raise ValueError("Unknown orbital_preference parameter.")
    else:
        state = ([1] + [0] * (numb_sites_per_particle - 1)) * int((system_size / numb_sites_per_particle))

    print("initial state = ", state)
    print("number of particles = ", numb_particles)
    print("number of lattice sites = ", len(state))
    if "Orbital" in model:
        print("number of orbital sites = ", 2 * len(state))

    return state


###################################################################
# define_iDMRG_model (defines the parameters for the iDMRG model) #
###################################################################


def define_iDMRG_model(model, ham_params):
    model_params = dict(conserve='N', bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle', verbose=1)
    basic_params = {k: ham_params[k] for k in ('t1', 'mu', 'V', 'Vtype', 'Vrange', 'n', 'nphi',  'LxMUC', 'Ly')}
    if 'phi' in ham_params:
        basic_params.update({'phi': ham_params['phi']})
    model_params.update(basic_params)

    if "Bos" in model:
        model_params.update(statistics='bosons', Nmax=1)
    else:  # "Fer"
        model_params.update(statistics='fermions')

    if 'HofSqu1' in model:
        M = HofSqu1Model(model_params)
    elif 'HofHex1' in model:
        M = HofHex1Model(model_params)
    elif 'HofHex1Hex5' in model:
        model_params.update(t2=model_params['t2'])
        M = HofHex1Hex5Model(model_params)
    else:  # "HofHex1Hex5Orbital"
        model_params.update(t2=model_params['t2'], t2dash=model_params['t2dash'], U=model_params['U'])
        M = HofHex1Hex5OrbitalModel(model_params)

    return M


#################################################################################
# my_iDMRG_pickle (defines an iDMRG engine or run iDMRG and determine pickling) #
#################################################################################


def my_iDMRG_pickle(program, model, chi_max, ham_params, use_pickle, make_pickle, run=True):

    # The run parameter specifies whether you are running iDMRG or defining an iDMRG engine. Defining an iDMRG engine
    # returns the engine, whereas running iDMRG returns [E, psi, M].

    pickle_file = None
    if use_pickle or make_pickle:
        if not run:
            pickle_stem = fp.file_name_stem("engine", model, chi_max)
        else:
            pickle_stem = fp.file_name_stem("E_psi_M", model, chi_max)
        pickle_leaf = fp.file_name_leaf("pickle", model, ham_params)
        # observables program needs to seek the output from ground_state
        os.makedirs(f"pickles/{program}/{model}/".replace("observables", "ground_state"), exist_ok=True)
        pickle_file = f"pickles/{program}/{model}/".replace("observables", "ground_state") + pickle_stem + pickle_leaf

    if use_pickle:
        with open(pickle_file, 'rb') as file1:
            if not run:
                engine = pickle.load(file1)
                shelve, sweep = False, 0  # shelved pickles not implemented for the flows
            else:
                [E, psi, M, shelve, sweep] = pickle.load(file1)
    else:
        shelve, sweep = False, 0

    if (program is not "observables" and shelve) or not use_pickle:
        if not use_pickle:
            engine = None
            (E, psi, M) = (None, None, None)
        if not run:
            engine = __my_iDMRG(model, chi_max, ham_params, shelve, sweep, run=False)
        else:
            (E, psi, M, shelve, sweep) = __my_iDMRG(model, chi_max, ham_params, shelve, sweep, run=True)
        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                if not run:
                    pickle.dump(engine, file2)
                else:
                    pickle.dump([E, psi, M, shelve, sweep], file2)

    if not run:
        return engine
    else:
        return E, psi, M


#####################################
# __my_iDMRG (implements the iDMRG) #
#####################################


def __my_iDMRG(model, chi_max, ham_params, shelve, sweep, run=True):

    M = define_iDMRG_model(model, ham_params)
    product_state = __get_product_state(model, ham_params)
    print(product_state)
    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        # 'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
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
        'verbose': 5,
        'N_sweeps_check': 10,
        'diag_method': 'lanczos',
        'max_hours': 14*24  # 2 weeks
    }

    if shelve:
        dmrg_params.update({'sweep_0': sweep})

    if not run:
        # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
        engine = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
        return engine
    else:
        info = dmrg.run(psi, M, dmrg_params)
        E = info['E']
        shelve = info['shelve']
        sweep = info['sweep_statistics']['sweep'][0]
        return E, psi, M, shelve, sweep


if __name__ == "__main__":

    __get_product_state(model="FerHofHex1Hex5Orbital", ham_params=dict(n=(1, 15), nphi=(1, 3), LxMUC=1, Ly=5),
                        filling_scale_factor=1, orbital_preference='polarizedy')
