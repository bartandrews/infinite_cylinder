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


##################################################
# get_product_state (calculates the initial psi) #
##################################################


def get_product_state(model, nn, nd, q, LxMUC, Ly, filling_scale_factor=1, orbital_preference=None):
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
    elif "Hex" in model:
        system_size = int(2 * q * LxMUC * Ly)
    else:
        raise ValueError("Unknown model for the get_product_state function.")
    numb_sites_per_particle = int(system_size / numb_particles)

    if "Orbital" in model:
        if "Bos" or "Fer" not in model:
            raise ValueError("Neither Bos nor Fer in model name.")
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


def define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                       nn, nd, p, q, LxMUC, Ly, phi=0):
    model_params = dict(conserve='N', t1=t1, mu=mu, V=V, Vtype=Vtype, Vrange=Vrange,
                        n=(int(nn), int(nd)), nphi=(int(p), int(q)),
                        LxMUC=LxMUC, Ly=Ly,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=phi)

    if "Bos" in model:
        model_params.update(statistics='bosons', Nmax=1)
    elif "Fer" in model:
        model_params.update(statistics='fermions')
    else:
        raise ValueError("Neither Bos nor Fer in model name.")

    if 'HofSqu1' in model:
        M = HofSqu1Model(model_params)
    elif 'HofHex1' in model:
        M = HofHex1Model(model_params)
    elif 'HofHex1Hex5' in model:
        model_params.update(t2=t2)
        M = HofHex1Hex5Model(model_params)
    elif 'HofHex1Hex5Orbital' in model:
        model_params.update(t2=t2, t2dash=t2dash, U=U)
        M = HofHex1Hex5OrbitalModel(model_params)
    else:
        raise ValueError("Unknown model for the define_iDMRG_model function.")

    if (V == 0 and Vrange != 0) or (V != 0 and Vrange == 0):
        raise ValueError("Cannot have zero interaction over a finite range, or a finite interaction over zero range.")

    return M


####################################################
# my_iDMRG (defines an iDMRG engine or runs iDMRG) #
####################################################


def my_iDMRG_pickle(flow, model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                    nn, nd, p, q, LxMUC, Ly,
                    use_pickle=False, make_pickle=False, phi=0, run=True):

    # The run parameter specifies whether you are running iDMRG or defining an iDMRG engine. Defining an iDMRG engine
    # returns the engine, whereas running iDMRG returns [E, psi, M].

    pickle_file = None
    if use_pickle or make_pickle:
        if not run:
            pickle_stem = fp.file_name_stem("engine", model, chi_max)
        else:
            pickle_stem = fp.file_name_stem("E_psi_M", model, chi_max)
        pickle_leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U}_mu_{mu}_V_{V}_{Vtype}_{Vrange}_" \
                      f"n_{nn}_{nd}_nphi_{p}_{q}_LxMUC_{LxMUC}_Ly_{Ly}_phi_{phi}.pkl"
        os.makedirs(f"pickles/{flow}/{model}/", exist_ok=True)
        pickle_file = f"pickles/{flow}/{model}/" + pickle_stem + pickle_leaf

    if use_pickle:
        with open(pickle_file, 'rb') as file1:
            if not run:
                engine = pickle.load(file1)
            else:
                [E, psi, M] = pickle.load(file1)
    else:
        engine = None
        (E, psi, M) = (None, None, None)
        if not run:
            engine = my_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange, nn, nd, p, q,
                              LxMUC, Ly, phi, run=False)
        else:
            (E, psi, M) = my_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange, nn, nd, p,
                                   q, LxMUC, Ly, phi, run=True)
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


def my_iDMRG(model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
             nn, nd, p, q,
             LxMUC, Ly, phi=0, run=True):

    M = define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                           nn, nd, p, q, LxMUC, Ly, phi)
    product_state = get_product_state(model, nn, nd, q, LxMUC, Ly)
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
    get_product_state(model="FerHofHex1Hex5Orbital", nn=1, nd=15, q=3, LxMUC=1, Ly=5,
                      filling_scale_factor=1, orbital_preference='polarizedy')
