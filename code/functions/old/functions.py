# --- python imports
import numpy as np
import sys
import pickle
# --- TeNPy imports
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
# from tenpy.algorithms.mps_sweeps import OneSiteH, TwoSiteH
# --- infinite_cylinder imports
from models.hofstadter import BosonicHofstadterModel, FermionicHofstadterModel
from models.hex_1 import BosonicHex1Model, FermionicHex1Model
from models.hex_1_hex_5 import BosonicHex1Hex5Model, FermionicHex1Hex5Model
from models.hex_1_hex_5_orbital import BosonicHex1Hex5OrbitalModel, FermionicHex1Hex5OrbitalModel


###############################################
# file_name_stem (creates the file name stem) #
###############################################


def file_name_stem(tool, model, chi_max):

    if model not in ['BosonicHofstadter', 'FermionicHofstadter',
                     'BosonicHex1', 'FermionicHex1',
                     'BosonicHex1Hex5', 'FermionicHex1Hex5',
                     'BosonicHex1Hex5Orbital', 'FermionicHex1Hex5Orbital']:
        sys.exit('Error: Unknown model.')

    stem = f"{tool}_{model}_chi_{chi_max}_"

    return stem


######################################################################
# Logger (writes stdout and stderr both to the screen and to a file) #
######################################################################


class Logger(object):
    def __init__(self, flow, model, leaf):
        self.terminal = sys.stdout or sys.stderr
        self.log = open(f"data/output/{flow}/{model}/output_{flow}_{model}_{leaf}", 'w', buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


######################################################################################
# print_LylB_headings (prints the Ly/lB ratio to file if there is a range of values) #
######################################################################################


def print_LylB_headings(model, Ly, ndvalue, nd_min, pvalue, qvalue, q_min, Ly_min):

    if "Hofstadter" in model:
        LylB = Ly * np.sqrt(2 * np.pi * (pvalue / qvalue))
    elif "Hex" in model:
        LylB = Ly * np.sqrt((4 * np.pi * (pvalue / qvalue)) / np.sqrt(3))
    else:
        sys.exit("Error: Unknown model for the print_LylB_headings function.")

    if ndvalue == nd_min and qvalue == q_min and Ly == Ly_min:
        data_line = f"LylB={LylB:.15f}\n"
    else:
        data_line = f"\n\nLylB={LylB:.15f}\n"

    return data_line


##################################################
# get_product_state (calculates the initial psi) #
##################################################


def get_product_state(model, nnvalue, ndvalue, qvalue, Lx_MUC, Ly, filling_scale_factor=1, orbital_preference=None):

    if "Orbital" in model:
        numb_particles = 2 * int(qvalue) * int(Lx_MUC) * int(Ly) * int(filling_scale_factor) * int(nnvalue) / int(
            ndvalue)
    else:
        numb_particles = int(qvalue) * int(Lx_MUC) * int(Ly) * int(filling_scale_factor) * int(nnvalue) / int(
            ndvalue)

    if not numb_particles.is_integer():
        sys.exit("Error: Cannot fit an integer number of particles into lattice geometry.")
    else:
        numb_particles = int(numb_particles)

    if "Hofstadter" in model:
        system_size = int(qvalue * Lx_MUC * Ly)
    elif "Hex" in model:
        system_size = int(2 * qvalue * Lx_MUC * Ly)
    numb_sites_per_particle = int(system_size / numb_particles)

    if "Orbital" in model:
        if orbital_preference in ['polarizedx', None]:
            if "Bosonic" in model:
                lattice_site = ['1_x 0_y']
            elif "Fermionic" in model:
                lattice_site = ['full_x empty_y']
            state = (lattice_site + [0] * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == 'polarizedy':
            if "Bosonic" in model:
                lattice_site = ['0_x 1_y']
            elif "Fermionic" in model:
                lattice_site = ['empty_x full_y']
            state = (lattice_site + [0] * (numb_sites_per_particle - 1)) * int(
                (system_size / numb_sites_per_particle))
        elif orbital_preference == 'unpolarized':
            if "Bosonic" in model:
                lattice_site_1 = ['1_x 0_y']
                lattice_site_2 = ['0_x 1_y']
            elif "Fermionic" in model:
                lattice_site_1 = ['full_x empty_y']
                lattice_site_2 = ['empty_x full_y']
            state = (lattice_site_1 + [0] * (numb_sites_per_particle - 1)
                     + lattice_site_2 + [0] * (numb_sites_per_particle - 1)) \
                    * int((system_size / (2 * numb_sites_per_particle)))
        elif orbital_preference == 'filled':
            if "Bosonic" in model:
                lattice_site = ['1_x 1_y']
            elif "Fermionic" in model:
                lattice_site = ['full_x full_y']
            state = (lattice_site + [0] * (2 * numb_sites_per_particle - 1)) \
                    * int((system_size / (2 * numb_sites_per_particle)))
        else:
            sys.exit("Error: Unknown orbital_preference parameter.")
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

def define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext=0):

    model_params = dict(conserve='N', t1=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)),
                        Lx_MUC=Lx_MUC, Ly=Ly,  # system params
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',  # MPS params
                        verbose=1, phi_ext=phi_ext)  # utility

    if "Bosonic" in model:
        model_params.update(Nmax=1)
    elif "Fermionic" in model:
        model_params.update(V=V)

    if model == 'BosonicHofstadter':
        M = BosonicHofstadterModel(model_params)
    elif model == 'FermionicHofstadter':
        M = FermionicHofstadterModel(model_params)
    elif model == 'BosonicHex1':
        M = BosonicHex1Model(model_params)
    elif model == 'FermionicHex1':
        M = FermionicHex1Model(model_params)
    elif model == 'BosonicHex1Hex5':
        model_params.update(t2=t2)
        M = BosonicHex1Hex5Model(model_params)
    elif model == 'FermionicHex1Hex5':
        model_params.update(t2=t2)
        M = FermionicHex1Hex5Model(model_params)
    elif model == 'BosonicHex1Hex5Orbital':
        model_params.update(t2=t2)
        model_params.update(t2dash=t2dash)
        model_params.update(U=U)
        M = BosonicHex1Hex5OrbitalModel(model_params)
    elif model == 'FermionicHex1Hex5Orbital':
        model_params.update(t2=t2)
        model_params.update(t2dash=t2dash)
        model_params.update(U=U)
        M = FermionicHex1Hex5OrbitalModel(model_params)

    return M


# #######################################################
# # define iDMRG (used when we want to reuse the state) #
# #######################################################
#
#
# def define_iDMRG_engine_pickle(flow, model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, use_pickle=False, make_pickle=False, phi_ext=0):
#
#     if use_pickle or make_pickle:
#         pickle_stem = file_name_stem("engine", model, chi_max)
#         pickle_leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_n_%s_%s_nphi_%s_%s_Lx_MUC_%s_Ly_%s_phi_%s.pkl" % (t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext))
#         pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf
#
#     if use_pickle:
#         with open(pickle_file, 'rb') as file1:
#             engine = pickle.load(file1)
#     else:
#         engine = define_iDMRG_engine(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
#         if make_pickle:
#             with open(pickle_file, 'wb') as file2:
#                 pickle.dump(engine, file2)
#
#     return engine
#
#
# def define_iDMRG_engine(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext=0):
#
#     M = define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
#     product_state = select_initial_psi(model, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly)
#     print(product_state)
#     psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)
#
#     dmrg_params = {
#         'mixer': True,  # setting this to True helps to escape local minima
#         # 'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
#         'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
#         'trunc_params': {
#             # 'chi_max': chi_max,
#             'svd_min': 1.e-10
#         },
#         # 'lanczos_params': {
#         #     'reortho': True,
#         #     'N_cache': 40
#         # },
#         'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},  # may need to change first value ~39
#         'max_E_err': 1.e-6,
#         'max_S_err': 1.e-6,
#         # 'norm_tol': 1.e-6,
#         # 'norm_tol_iter': 1000,
#         'max_sweeps': 1000,
#         'verbose': 5,
#         'N_sweeps_check': 10,
#         'diag_method': 'lanczos'
#     }
#
#     # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
#     engine = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
#
#     return engine
#
#
# ##########################################################
# # run iDMRG (used when we want to recalculate the state) #
# ##########################################################
#
#
# def run_iDMRG_pickle(flow, model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, use_pickle=False, make_pickle=False, phi_ext=0):
#
#     if use_pickle or make_pickle:
#         pickle_stem = file_name_stem("E_psi_M", model, chi_max)
#         pickle_leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_n_%s_%s_nphi_%s_%s_Lx_MUC_%s_Ly_%s_phi_%s.pkl" % (t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext))
#         pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf
#
#     if use_pickle:
#         with open(pickle_file, 'rb') as file1:
#             [E, psi, M] = pickle.load(file1)
#     else:
#         (E, psi, M) = run_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
#         if make_pickle:
#             with open(pickle_file, 'wb') as file2:
#                 pickle.dump([E, psi, M], file2)
#
#     return E, psi, M
#
#
# def run_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext=0):
#
#     M = define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
#     product_state = select_initial_psi(model, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly)
#     print(product_state)
#     psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)
#
#     dmrg_params = {
#         'mixer': True,  # setting this to True helps to escape local minima
#         'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
#         'trunc_params': {
#             # 'chi_max': chi_max,
#             'svd_min': 1.e-10
#         },
#         # 'lanczos_params': {
#         #     'reortho': True
#         #     'N_cache': 40
#         # },
#         'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
#         'max_E_err': 1.e-6,
#         'max_S_err': 1.e-6,
#         'max_sweeps': 1000,
#         'verbose': 1,
#         'N_sweeps_check': 10,
#         # 'diag_method': 'lanczos'
#     }
#
#     info = dmrg.run(psi, M, dmrg_params)
#     E = info['E']
#
#     return E, psi, M


####################################################
# my_iDMRG (defines an iDMRG engine or runs iDMRG) #
####################################################


def my_iDMRG_pickle(flow, model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly,
                    use_pickle=False, make_pickle=False, phi_ext=0, run=True):

    # The run parameter specifies whether you are running iDMRG or defining an iDMRG engine. Defining an iDMRG engine
    # returns the engine, whereas running iDMRG returns [E, psi, M].

    if use_pickle or make_pickle:
        if not run:
            pickle_stem = file_name_stem("engine", model, chi_max)
        else:
            pickle_stem = file_name_stem("E_psi_M", model, chi_max)
        pickle_leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_n_%s_%s_nphi_%s_%s_Lx_MUC_%s_Ly_%s_phi_%s.pkl" % (t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:
        with open(pickle_file, 'rb') as file1:
            if run == False:
                engine = pickle.load(file1)
            else:
                [E, psi, M] = pickle.load(file1)
    else:
        if not run:
            engine = my_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue,
                                         Lx_MUC, Ly, phi_ext, run=False)
        else:
            (E, psi, M) = my_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC,
                                    Ly, phi_ext, run=True)
        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                if not run:
                    pickle.dump(engine, file2)
                else:
                    pickle.dump([E, psi, M], file2)

    if not run:
        return engine
    else:
        return E, psi, M


def my_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext=0, run=True):

    M = define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
    product_state = get_product_state(model, nnvalue, ndvalue, qvalue, Lx_MUC, Ly)
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
        #     'reortho': True,
        #     'N_cache': 40
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
        'diag_method': 'lanczos'
    }

    if not run:
        # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
        engine = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
        return engine
    else:
        info = dmrg.run(psi, M, dmrg_params)
        E = info['E']
        return E, psi, M


if __name__ == "__main__":

    get_product_state(model="FermionicHex1Hex5Orbital", nnvalue=1, ndvalue=15, qvalue=3, Lx_MUC=1, Ly=5,
                      filling_scale_factor=1, orbital_preference='polarizedy')
