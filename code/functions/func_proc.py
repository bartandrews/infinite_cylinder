# --- python imports
import sys
import os
import inspect  # for main
import pkgutil  # for main
# --- infinite_cylinder imports
from models.hofstadter.hofstadter import HofstadterModel  # for main
import models.hofstadter as hofstadter  # for main


#######################################################
# check_input_params (check for errors in user input) #
#######################################################


def check_input_params(program, threads, model, chi_max, ham_params, use_pickle=False, make_pickle=False):

    if not isinstance(threads, int):
        raise TypeError("threads needs to be an integer.")
    if threads < 0:
        raise ValueError("threads needs to be positive.")

    if model not in ['BosHofSqu1', 'FerHofSqu1',
                     'BosHofHex1', 'FerHofHex1',
                     'BosHofHex1Hex5', 'FerHofHex1Hex5',
                     'BosHofHex1Hex5Orbital', 'FerHofHex1Hex5Orbital']:
        raise ValueError('Unknown model for the check_input_params function.')

    if not isinstance(chi_max, int):
        raise TypeError("chi_max needs to be an integer.")
    if chi_max < 0:
        raise ValueError("chi_max needs to be positive.")

    if not isinstance(use_pickle, bool) or not isinstance(make_pickle, bool):
        raise TypeError("use_pickle and make_pickle need to be boolean.")

    flow_variables = ['phi', 'U', 'V', 'kappa']
    for i in flow_variables:
        if program is f'{i}_flow':
            if f'{i}_min' or f'{i}_max' or f'{i}_samp' not in ham_params:
                raise KeyError(f"{i}_min, {i}_max, {i}_samp need to be specified for {i}_flow.")
            if ham_params[f'{i}_min'] > ham_params[f'{i}_max']:
                raise ValueError(f"{i}_max has to be greater than {i}_min.")
            if isinstance(ham_params[f't{i}_min'], str) or isinstance(ham_params[f't{i}_max'], str) \
                    or isinstance(ham_params[f't{i}_samp'], str):
                raise TypeError(f"{i}_min, {i}_max, {i}_samp need to be numbers.")

    for i in range(1, 11, 1):  # search up to 10th-NN hoppings for both t and tdash
        if f"{i}" in model:
            if isinstance(ham_params[f't{i}'], str):
                raise TypeError("All t hoppings need to be numbers.")
            if f"t{i}dash" in ham_params:
                if isinstance(ham_params[f't{i}dash'], str):
                    raise TypeError("All tdash hoppings need to be numbers.")

    if 'mu' in ham_params:
        if isinstance(ham_params['mu'], str):
            raise TypeError("mu needs to be a number.")

    if 'V' in ham_params:
        if isinstance(ham_params['V'], str):
            raise TypeError("V needs to be a number.")

    if 'Vtype' in ham_params:
        if ham_params['Vtype'] not in ['Coulomb', 'Yukawa']:
            raise ValueError('Unknown Vtype for the check_input_params function.')

    if 'Vrange' in ham_params:
        if ham_params['Vrange'] not in range(11):
            raise ValueError("Vrange needs to be an integer in [0, 10].")

    if 'V' in ham_params:
        if (ham_params['V'] == 0 and ham_params['Vrange'] != 0) \
                or (ham_params['V'] != 0 and ham_params['Vrange'] == 0):
            raise ValueError("Cannot have zero interaction over a finite range, "
                             "or a finite interaction over zero range.")

    if 'n' in ham_params:
        if not isinstance(ham_params['n'][0], int) or not isinstance(ham_params['n'][1], int):
            raise TypeError("n needs to have integer entries.")
        if ham_params['n'][0] < 0 or ham_params['n'][1] < 0:
            raise ValueError("n needs to have positive entries.")

    if 'nphi' in ham_params:
        if not isinstance(ham_params['nphi'][0], int) or not isinstance(ham_params['nphi'][1], int):
            raise TypeError("nphi needs to have integer entries.")
        if ham_params['nphi'][0] < 0 or ham_params['nphi'][1] < 0:
            raise ValueError("nphi needs to have positive entries.")

    if 'LxMUC' and 'Ly' in ham_params:
        if not isinstance(ham_params['LxMUC'], int) or not isinstance(ham_params['Ly'], int):
            raise TypeError("LxMUX and Ly need to be integers.")
        if ham_params['LxMUC'] < 0 or ham_params['Ly'] < 0:
            raise ValueError("LxMUX and Ly need to be positive.")

    if 'phi' in ham_params:
        if isinstance(ham_params['phi'], str):
            raise TypeError("phi needs to be a number.")

    if 'tag' in ham_params:
        if not isinstance(ham_params['tag'], str):
            raise TypeError("tag needs to be a string.")


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

    if program is "kappa_flow":
        kappa = f"kappa_{ham_params['kappa_min']:g}_{ham_params['kappa_max']:g}_{ham_params['kappa_samp']}_"
    else:
        kappa = ""

    if program is "U_flow":
        U = f"U_{ham_params['U_min']:g}_{ham_params['U_max']:g}_{ham_params['U_samp']}_"
    else:
        U = f"U_{ham_params['U']:g}_" if ham_params['U'] != 0 else ""

    mu = f"mu_{ham_params['mu']:g}_" if ham_params['mu'] != 0 else ""

    if program is "V_flow":
        V = f"V_{ham_params['V_min']:g}_{ham_params['V_max']:g}_{ham_params['V_samp']}_{ham_params['Vtype']}_{ham_params['Vrange']}_"
    else:
        V = f"V_{ham_params['V']:g}_{ham_params['Vtype']}_{ham_params['Vrange']}_" if ham_params['V'] != 0 else ""

    nu = f"n_{ham_params['n'][0]}_{ham_params['n'][1]}_nphi_{ham_params['nphi'][0]}_{ham_params['nphi'][1]}_"
    L = f"LxMUC_{ham_params['LxMUC']}_Ly_{ham_params['Ly']}"

    if program is "phi_flow":
        phi = f"_phi_{ham_params['phi_min']:g}_{ham_params['phi_max']:g}_{ham_params['phi_samp']}"
    else:
        if "phi" in ham_params:
            phi = f"phi_{ham_params['phi']:g}" if ham_params['phi'] != 0 else ""
        else:
            phi = ""

    ext = ".dat" if program is not "pickle" else ".pkl"

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
