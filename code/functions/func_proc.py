# --- python imports
import sys
import os
import ntpath
import inspect  # for main
import pkgutil  # for main
# --- infinite_cylinder imports
from models.hofstadter.hofstadter import HofstadterModel  # for main
import models.hofstadter as hofstadter  # for main


#######################################################################
# process_pickle_file_name (extract parameters from pickle file name) #
#######################################################################


def process_pickle_file_name(filepath):

    if ".pkl" not in filepath:
        raise ValueError("pickle file needs to parsed as first argument.")
    if "E_psi_M" not in filepath:
        raise ValueError("pickle file of the ground state needs to be targeted.")

    print("filepath = ", filepath)
    pickle = __file_name(filepath)
    print("pickle = ", pickle)
    debased_pickle = str(pickle.replace("E_psi_M_", "").split(".pkl", 1)[0])
    print("debased_pickle = ", debased_pickle)
    debased_pickle_entries = debased_pickle.split('_')
    model = debased_pickle_entries[0]
    print("model = ", model)
    chi_max = int(debased_pickle_entries[2])
    print("chi_max = ", chi_max)
    LxMUC = int(debased_pickle_entries[debased_pickle_entries.index("LxMUC")+1])
    print("LxMUC = ", LxMUC)
    Ly = int(debased_pickle_entries[debased_pickle_entries.index("Ly") + 1])
    print("Ly = ", Ly)
    leaf_entries = debased_pickle_entries[3:]
    leaf = '_'.join(leaf_entries) + ".dat"
    print("leaf = ", leaf)

    return model, chi_max, leaf, LxMUC, Ly


#############################################
# __file_name (extract file name from path) #
#############################################


def __file_name(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


###############################################
# file_name_leaf (creates the file name leaf) #
###############################################


def file_name_leaf(program, model, ham_params):

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

    nu = f"n_{ham_params['n'][0]}_{ham_params['n'][1]}_nphi_{ham_params['nphi'][0]}_{ham_params['nphi'][1]}_"
    L = f"LxMUC_{ham_params['LxMUC']}_Ly_{ham_params['Ly']}"

    if program == "phi_flow":
        phi = f"_phi_{ham_params['phi_min']:g}_{ham_params['phi_max']:g}_{ham_params['phi_samp']}"
    else:
        if "phi" in ham_params:
            phi = f"phi_{ham_params['phi']:g}" if ham_params['phi'] != 0 else ""
        else:
            phi = ""

    ext = ".dat" if program != "pickle" else ".pkl"

    leaf = f"{t}{kappa}{U}{mu}{V}{nu}{L}{phi}{ext}{ham_params['tag']}"

    return leaf


###############################################
# file_name_stem (creates the file name stem) #
###############################################


def file_name_stem(tool, model, chi_max):
    return f"{tool}_{model}_chi_{chi_max}_".replace(" ", "_")


###################################################################
# prepare_output_files (creates the output directories and files) #
###################################################################


def prepare_output_files(tools, model, chi_max, leaf, chiK_max=0):

    stem, file, data = [dict()]*3
    for tool in tools:
        stem.update({tool: file_name_stem(tool, model, chi_max)})
        os.makedirs(f"data/{tool}/{model}/", exist_ok=True)
        if tool == 'ent_spec_mom':
            file.update({tool: f"data/{tool}/{model}/" + stem[tool] + f"chiK_{chiK_max}_" + leaf})
        else:
            file.update({tool: f"data/{tool}/{model}/" + stem[tool] + leaf})
        open(file[tool], "w")
        data[tool] = open(file[tool], "a", buffering=1)

    return data


######################################################################
# Logger (writes stdout and stderr both to the screen and to a file) #
######################################################################


class Logger(object):
    def __init__(self, program, model, chi_max, leaf):
        self.terminal = sys.stdout or sys.stderr
        os.makedirs(f"logs/{program}/{model}/", exist_ok=True)
        self.log = open(f"logs/{program}/{model}/log_{program}_{model}_chi_{chi_max}_{leaf}", 'w', buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


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
