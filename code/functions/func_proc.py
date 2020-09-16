# --- python imports
import sys
import os
import inspect  # for main
import pkgutil  # for main
import fnmatch
import binascii
# --- infinite_cylinder imports
from models.hofstadter.hofstadter import HofstadterModel  # for main
import models.hofstadter as hofstadter  # for main


###############################################
# file_name_stem (creates the file name stem) #
###############################################


def file_name_stem(tool, model, chi_max):
    return f"{tool}_{model}_chi_{chi_max}_".replace(" ", "_")


###############################################
# file_name_leaf (creates the file name leaf) #
###############################################


def file_name_leaf(program, model, ham_params):

    if "Nmax" in ham_params:
        Nmax = f"Nmax_{ham_params['Nmax']}_" if ham_params['Nmax'] != 1 else ""
    else:
        Nmax = ""

    C = f"C_{ham_params['C']}_" if "HalSquCN" in model else ""

    if program == "ratio_flow":
        ratio = f"r_{ham_params['r_min']:g}_{ham_params['r_max']:g}_{ham_params['r_samp']}_"
    else:
        if "r" in ham_params:
            ratio = f"r_{ham_params['r']:g}" if ham_params['r'] is not None else ""
        else:
            ratio = ""

    t = ""
    for i in range(1, 11, 1):  # search up to 10th-NN hoppings for both t and tdash
        if f"{i}" in model:
            t += f"t{i}_{ham_params[f't{i}']:g}_" if ham_params[f't{i}'] != 0 else ""
            if f"t{i}dash" in ham_params:
                t += f"t{i}dash_{ham_params[f't{i}dash']:g}_" if ham_params[f't{i}dash'] != 0 else ""

    if program == "kappa_flow":
        kappa = f"kappa_{ham_params['kappa_min']:g}_{ham_params['kappa_max']:g}_{ham_params['kappa_samp']}_"
    else:
        kappa = ""

    if program == "U_flow":
        U = f"U_{ham_params['U_min']:g}_{ham_params['U_max']:g}_{ham_params['U_samp']}_"
    else:
        U = f"U_{ham_params['U']:g}_" if ham_params['U'] != 0 else ""

    mu = f"mu_{ham_params['mu']:g}_" if ham_params['mu'] != 0 else ""

    if program == "V_flow":
        V = f"V_{ham_params['V_min']:g}_{ham_params['V_max']:g}_{ham_params['V_samp']}_{ham_params['Vtype']}_{ham_params['Vrange']}_"
    else:
        V = f"V_{ham_params['V']:g}_{ham_params['Vtype']}_{ham_params['Vrange']}_" if ham_params['V'] != 0 else ""

    if "Hof" in model:
        nu = f"n_{ham_params['n'][0]}_{ham_params['n'][1]}_nphi_{ham_params['nphi'][0]}_{ham_params['nphi'][1]}_"
    else:
        nu = f"n_{ham_params['n'][0]}_{ham_params['n'][1]}_"
    L = f"LxMUC_{ham_params['LxMUC']}_Ly_{ham_params['Ly']}"

    if program == "phi_flow":
        phi = f"_phi_{ham_params['phi_min']:g}_{ham_params['phi_max']:g}_{ham_params['phi_samp']}"
    else:
        if "phi" in ham_params:
            phi = f"phi_{ham_params['phi']:g}" if ham_params['phi'] != 0 else ""
        else:
            phi = ""

    custom = "_custom" if ham_params['custom'] else ""

    if program == "h5":
        ext = ".h5"
    elif program == "pkl":
        ext = ".pkl"
    else:
        ext = ".dat"

    leaf = f"{Nmax}{C}{ratio}{t}{kappa}{U}{mu}{V}{nu}{L}{phi}{custom}{ext}{ham_params['tag']}"

    return leaf


######################################################################
# Logger (writes stdout and stderr both to the screen and to a file) #
######################################################################


class Logger(object):
    def __init__(self, program, path, model, chi_max, leaf):
        self.terminal = sys.stdout or sys.stderr
        os.makedirs(os.path.join(path, "logs", f"{program}", f"{model}", ""), exist_ok=True)
        self.log = open(os.path.join(path, "logs", f"{program}", f"{model}",
                                     f"log_{program}_{model}_chi_{chi_max}_{leaf}"), 'w', buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


###################################################################
# prepare_output_files (creates the output directories and files) #
###################################################################


def prepare_output_files(tools, path, model, chi_max, leaf, chiK_max=0):

    stem, file, data = [dict()]*3
    for tool in tools:
        stem.update({tool: file_name_stem(tool, model, chi_max)})
        os.makedirs(os.path.join(path, "data", f"{tool}", f"{model}", ""), exist_ok=True)
        if tool == 'ent_spec_mom':
            file.update({tool: os.path.join(path, "data", f"{tool}", f"{model}",
                                            stem[tool] + f"chiK_{chiK_max}_" + leaf)})
        else:
            file.update({tool: os.path.join(path, "data", f"{tool}", f"{model}", stem[tool] + leaf)})
        open(file[tool], "w")
        data[tool] = open(file[tool], "a", buffering=1)

    return data

###########################################
# is_gz_file (check if a file is gzipped) #
###########################################


def is_gz_file(my_file_path):
    with open(my_file_path, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


#############################################################
# process_pickle (extract parameters from pickle file name) #
#############################################################


def process_pickle(my_file):

    extra_dof_list = ["Orbital", "Spin"]

    if "state" in my_file:
        prefix = "state_"
    elif "E_psi_M" in my_file:  # backward compatibility
        prefix = "E_psi_M_"
    else:
        raise ValueError("pickle file of the ground state needs to be targeted.")

    if ".h5" in my_file:
        ext = ".h5"
    elif ".pkl" in my_file:  # backward compatibility
        ext = ".pkl"
    else:
        raise ValueError("pickle file needs to parsed as first argument.")

    print("pickle = ", my_file)
    debased_pickle = str(my_file.replace(prefix, "").split(ext, 1)[0])
    print("debased_pickle = ", debased_pickle)
    debased_pickle_entries = debased_pickle.split('_')
    model = debased_pickle_entries[0]
    print("model = ", model)
    chi_max = int(debased_pickle_entries[2])
    print("chi_max = ", chi_max)
    leaf_entries = debased_pickle_entries[3:]
    leaf = '_'.join(leaf_entries) + ".dat"
    print("leaf = ", leaf)
    extra_dof_flag = True if any(dof in model for dof in extra_dof_list) else False
    print("extra_dof_flag = ", extra_dof_flag)

    return model, chi_max, leaf, extra_dof_flag


###############################################################################################
# largest_chi_pickle (for a given configuration, return the pickle path with the largest chi) #
###############################################################################################


def largest_chi_pickle(my_file_path, chi_max):

    (pickle_dir, pickle_file) = os.path.split(my_file_path)

    # get the pkl file list for a given configuration at various chi
    complete_list = os.listdir(path=pickle_dir)

    file_split = pickle_file.split('_')
    file_split[file_split.index("chi") + 1] = "*"
    file_join = '_'.join(file_split)
    file_join_2 = file_join.replace('.pkl', '.*').replace('.h5', '.*')
    file_general = file_join_2  # tagged files are not processed

    pkl_file_list = []
    for i, val in enumerate(complete_list):  # iterate through all files in the pickle directory
        if fnmatch.fnmatch(val, file_general):
            pkl_file_list.append(val)

    # identify the file with largest chi for the pkl file list for a given configuration
    largest_index = 0
    largest_chi = 0
    for i, val in enumerate(pkl_file_list):  # iterate through all files in the pickle directory
        pkl1 = val.split('_')
        if int(pkl1[pkl1.index('chi') + 1]) > largest_chi:
            largest_index = i
            largest_chi = int(pkl1[pkl1.index('chi') + 1])
    largest_pkl_file = pkl_file_list[largest_index]

    if chi_max == largest_chi:
        print("Desired chi is equal to largest available chi.")
        target_pickle_path = my_file_path
    elif chi_max > largest_chi:
        print("Desired chi is larger than largest available chi.")
        target_pickle_path = os.path.join(pickle_dir, largest_pkl_file)
    else:  # chi_max < largest_chi
        raise ValueError("Desired chi is smaller than largest available chi. Pickle cannot be used.")

    return target_pickle_path


if __name__ == '__main__':

    # lists the submodules but not the subclasses
    for importer, modname, ispkg in pkgutil.iter_modules(hofstadter.__path__):
        print("Found submodule %s (is a package: %s)" % (modname, ispkg))
        print(inspect.getmembers(sys.modules[__name__], inspect.isclass))

    models = []
    for m in HofstadterModel.__subclasses__():  # only shows the subclasses of imported modules
        print(m)
        models.append(m.__name__)
    print(models)
