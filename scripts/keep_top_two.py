# keep_top_two.py --- delete all files other than those with the largest two values of chi for each configuration (execute in directory)
#
# Conditions:
# - tagged files are ignored
# - this script assumes the same number of entries in each file name

# --- python imports
import os
import re
from itertools import groupby
from fractions import Fraction as Frac
from operator import itemgetter
import ast


###################################
# key function for the model sort #
###################################


def model_keyf(a):
    start_string = re.split('_chi_(\d+)_', a)[0]
    start_string_1 = start_string.replace("state_", "")
    start_string_2 = start_string_1.replace("E_psi_M_", "")
    return start_string_2


####################################
# key function for the system sort #
####################################


# system_keyf = lambda a: re.split('_chi_(\d+)_', a)[2]
def system_keyf(a):
    return re.split('_chi_(\d+)_', a)[2]


##############################
# auto convert types in list #
##############################


def tryeval(value):
    try:
        value = ast.literal_eval(value)
    except ValueError:
        pass
    return value


###############################
# decompose the dat filenames #
###############################


def sort_list(mylist):

    # decompose list
    decomposed_list = []
    for _i, _val in enumerate(mylist):
        if _val.startswith('log_observables_'):
            prefix = "log_observables_"
            suffix = ".dat"
        elif _val.startswith('log_ground_state_'):
            prefix = "log_ground_state_"
            suffix = ".dat"
        elif _val.startswith('state_'):
            prefix = "state_"
            suffix = ".pkl"
        elif _val.startswith('E_psi_M_'):
            prefix = "E_psi_M_"
            suffix = ".pkl"
        else:
            raise ValueError("Unknown prefix for the file: ", _val)
        decomposed_list.append(str(_val.replace(prefix, "").split(suffix, 1)[0]).split('_'))

    # add nu entries to the end of each decomposed configuration
    for _i, _val in enumerate(mylist):
        _nn = int(decomposed_list[_i][decomposed_list[_i].index("n") + 1])
        _nd = int(decomposed_list[_i][decomposed_list[_i].index("n") + 2])
        _p = int(decomposed_list[_i][decomposed_list[_i].index("nphi") + 1])
        _q = int(decomposed_list[_i][decomposed_list[_i].index("nphi") + 2])
        _nphi = _p / _q
        _nu = (_nn / _nd) / _nphi
        _frac_nu = Frac(str(_nu)).limit_denominator(100)
        decomposed_list[_i] += ['nu', _frac_nu.numerator, _frac_nu.denominator]

    # auto convert the types of entries in the list
    for _i, _val in enumerate(mylist):
        decomposed_list[_i] = [tryeval(x) for x in decomposed_list[_i]]

    # sort list
    model_index = 0
    r_index = decomposed_list[0].index("nu") + 1
    s_index = decomposed_list[0].index("nu") + 2
    p_index = decomposed_list[0].index("nphi") + 1
    q_index = decomposed_list[0].index("nphi") + 2
    Ly_index = decomposed_list[0].index("Ly") + 1
    chi_index = decomposed_list[0].index("chi") + 1

    sorted_list = sorted(decomposed_list, key=itemgetter(model_index, r_index, s_index,
                                                         p_index, q_index, Ly_index, chi_index))

    # convert the types of entries in the list back into strings
    for _i, _val in enumerate(mylist):
        sorted_list[_i] = [str(x) for x in sorted_list[_i]]

    # remove nu entries and recompose list
    recomposed_list = []
    for _i, _val in enumerate(sorted_list):
        sorted_list[_i] = sorted_list[_i][:len(sorted_list[_i]) - 3]
        if os.path.exists("log_observables_" + '_'.join(sorted_list[_i]) + ".dat"):
            recomposed_list.append("log_observables_" + '_'.join(sorted_list[_i]) + ".dat")
        elif os.path.exists("log_ground_state_" + '_'.join(sorted_list[_i]) + ".dat"):
            recomposed_list.append("log_ground_state_" + '_'.join(sorted_list[_i]) + ".dat")
        elif os.path.exists("state_" + '_'.join(sorted_list[_i]) + ".pkl"):
            recomposed_list.append("state_" + '_'.join(sorted_list[_i]) + ".pkl")
        elif os.path.exists("E_psi_M_" + '_'.join(sorted_list[_i]) + ".pkl"):
            recomposed_list.append("E_psi_M_" + '_'.join(sorted_list[_i]) + ".pkl")
        else:
            raise ValueError("Missing file corresponding to the configuration: ", '_'.join(sorted_list[_i]))

    return recomposed_list


if __name__ == '__main__':

    # get the complete directory list
    complete_dir_list = os.listdir()

    # get the (dat or pkl) file list
    file_list = []
    suffixes = [".dat", ".pkl"]
    for i, val in enumerate(complete_dir_list):
        if any(x in val for x in suffixes):
            if "Nmax" in val:
                continue
            else:
                file_list.append(val)

    # sort the file list
    sorted_file_list = sort_list(file_list)

    # group by string before the "_chi_*_" omitting "state_"/"E_psi_M_" prefix (by model)
    model_grouped_list = [list(i) for j, i in groupby(sorted_file_list, key=model_keyf)]

    # subgroup by string after the "_chi_*_" (by system)
    system_grouped_list = []
    for i in range(len(model_grouped_list)):
        system_grouped_list.append(
            [list(i) for j, i in groupby(model_grouped_list[i], key=system_keyf)])

    # construct a list of the pickles that are actually used
    used_list = []
    for i in range(len(system_grouped_list)):  # loop over models
        for j in range(len(system_grouped_list[i])):  # loop over systems
            if len(system_grouped_list[i][j]) >= 2:
                used_list.append(system_grouped_list[i][j][-2:])
            elif len(system_grouped_list[i][j]) == 1:
                used_list.append(system_grouped_list[i][j][-1:])
            else:
                continue

    used_list = sum(used_list, [])  # flatten the used list

    list_difference = [item for item in file_list if item not in used_list]

    for val in list_difference:
        print(val)

    if input("Are you sure you want to remove the above files? (y/n): ") == "y":
        for val in list_difference:
            os.remove(val)
    else:
        exit()
