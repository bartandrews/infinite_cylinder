# --- python imports
import os
import pickle
# --- TeNPy imports
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
# from tenpy.algorithms.mps_sweeps import OneSiteH, TwoSiteH
# --- infinite_cylinder imports
import functions.func_proc as fp
from models.haldane.haldane import HalModel
from models.hofstadter.squ_1 import HofSqu1Model
from models.hofstadter.hex_1 import HofHex1Model
from models.hofstadter.hex_1_hex_5 import HofHex1Hex5Model
from models.hofstadter.hex_1_hex_5_orbital import HofHex1Hex5OrbitalModel
from models.old.magnetic_lattice.hex_1_hex_5_orbital import FermionicHex1Hex5OrbitalModel  # old code


####################################################
# __get_custom_state (define a custom initial psi) #
####################################################

def __get_custom_state():

    state = [1, 0]*3

    return state


####################################################
# __get_product_state (calculates the initial psi) #
####################################################


def __get_product_state(model, ham_params, filling_scale_factor=1, orbital_preference=None):

    nn, nd, LxMUC, Ly = ham_params['n'][0], ham_params['n'][1], ham_params['LxMUC'], ham_params['Ly']
    q = ham_params['nphi'][1] if "Hof" in model else 1

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
        empty_site = [0]
        state = ([1] + empty_site * (numb_sites_per_particle - 1)) * int((system_size / numb_sites_per_particle))

    state += (system_size-len(state)) * empty_site  # top up remainder with empty sites, if necessary

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
    basic_params = {k: ham_params[k] for k in ('Nmax', 't1', 'mu', 'V', 'Vtype', 'Vrange', 'n', 'nphi',  'LxMUC', 'Ly')}
    if 'phi' in ham_params:
        basic_params.update({'phi': ham_params['phi']})
    model_params.update(basic_params)

    if "Bos" in model:
        model_params.update(statistics='bosons')
    else:  # "Fer"
        model_params.update(statistics='fermions')

    if model.endswith("Hal"):
        del model_params['Vrange'], model_params['Vtype'], model_params['n'], model_params['nphi']
        M = HalModel(model_params)
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


def my_iDMRG_pickle(program, path, model, chi_max, ham_params, use_pickle, make_pickle, run=True):

    # The run parameter specifies whether you are running iDMRG or defining an iDMRG engine. Defining an iDMRG engine
    # returns the engine, whereas running iDMRG returns [E, psi, M].

    pickle_path = None
    if use_pickle or make_pickle:  # determine the given pickle path
        if run:
            pickle_stem = fp.file_name_stem("E_psi_M", model, chi_max)
        else:
            pickle_stem = fp.file_name_stem("engine", model, chi_max)
        pickle_leaf = fp.file_name_leaf("pickle", model, ham_params)
        os.makedirs(os.path.join(path, "pickles", f"{program}", f"{model}", ""), exist_ok=True)
        pickle_dir = os.path.join(path, "pickles", f"{program}", f"{model}")
        pickle_file = pickle_stem + pickle_leaf
        pickle_path = os.path.join(pickle_dir, pickle_file)

    improve_flag = False
    if use_pickle:  # get the information from the pickle
        target_pickle_path = fp.largest_chi_pickle(pickle_dir, pickle_file, pickle_path, chi_max)
        with open(target_pickle_path, 'rb') as file1:
            if run:  # sensible possibilities: resuming a shelved pickle or improving on a complete pickle
                print("Reading pickle: ", target_pickle_path)
                print("Writing pickle: ", pickle_path)
                [E, psi, M, shelve, sweep] = pickle.load(file1)
                if target_pickle_path != pickle_path:  # improving on a complete pickle
                    if shelve:
                        raise ValueError("Cannot improve on a shelved pickle. Complete the pickle first by running at the same chi.")
                    else:
                        improve_flag = True  # mark that we are improving on a complete pickle
            else:  # program == flow
                engine = pickle.load(file1)
                shelve, sweep = False, 0  # shelved pickles not implemented for the flows
    else:
        engine = None
        (E, psi, M, shelve, sweep) = (None, None, None, False, 0)

    if (program != "observables" and shelve) or not use_pickle or improve_flag:  # get the information by running iDMRG (either the information that we have is not complete, or we do not have the information, or we want to improve on the information)
        if run:
            if improve_flag:  # improving on a complete pickle
                (E, psi, M, shelve, sweep) = __my_iDMRG(model, chi_max, ham_params, shelve, sweep, run=True, initial_M=M, initial_psi=psi)
            else:
                (E, psi, M, shelve, sweep) = __my_iDMRG(model, chi_max, ham_params, shelve, sweep, run=True)
        else:
            engine = __my_iDMRG(model, chi_max, ham_params, shelve, sweep, run=False)
        if make_pickle:
            with open(pickle_path, 'wb') as file2:
                if not run:
                    pickle.dump(engine, file2)
                else:
                    pickle.dump([E, psi, M, shelve, sweep], file2)

    if run:
        return E, psi, M, shelve, sweep
    else:
        return engine


#####################################
# __my_iDMRG (implements the iDMRG) #
#####################################


def __my_iDMRG(model, chi_max, ham_params, shelve, sweep, run=True, initial_M=None, initial_psi=None):

    if (initial_M is not None) and (initial_psi is not None):  # improving on a complete pickle
        M, psi = initial_M, initial_psi
    else:
        M = define_iDMRG_model(model, ham_params)
        product_state = __get_product_state(model, ham_params) if not ham_params['custom'] else __get_custom_state()
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
        'verbose': 5,  # 5
        'N_sweeps_check': 10,
        'diag_method': 'default',
        'max_hours': 24*14  # 2 weeks
    }

    if model == "FerHofHex1Hex5OrbitalOld":  # old code
        del dmrg_params['lanczos_params'], dmrg_params['max_sweeps'], dmrg_params['diag_method'], \
            dmrg_params['max_hours']
        dmrg_params.update(mixer_params={'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30}, max_E_err=1.e-10,
                           verbose=1)

    if shelve:
        dmrg_params.update({'sweep_0': sweep})

    if (initial_M is not None) and (initial_psi is not None):  # improving on a complete pickle
        dmrg_params['mixer'] = False
        del dmrg_params['mixer_params']
        dmrg_params['trunc_params']['chi_max'] = chi_max
        del dmrg_params['chi_list']
        # dmrg_params.update({'sweep_0': sweep})

    if run:
        info = dmrg.run(psi, M, dmrg_params)
        E = info['E']
        shelve = info['shelve']
        sweep = info['sweep_statistics']['sweep'][-1]
        return E, psi, M, shelve, sweep
    else:
        # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
        engine = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
        return engine


if __name__ == "__main__":

    __get_product_state(model="FerHofSqu1", ham_params=dict(n=(4, 45), nphi=(4, 15), LxMUC=1, Ly=6),
                        filling_scale_factor=1, orbital_preference=None)
