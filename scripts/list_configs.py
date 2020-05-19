# list_configs.py --- generate a sorted list of all configurations in a directory (execute in directory)

import os
import re
from fractions import Fraction as Frac
from operator import itemgetter
import ast


###################################
# key function for the model sort #
###################################


# model_keyf = lambda a: re.split('_chi_(\d+)_', a)[0]
def model_keyf(a):
    return re.split('_chi_(\d+)_', a)[0]


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
            recomp_flag = "log_observables_"
            recomp_suffix = ".dat"
            decomposed_list.append(str(_val.replace(recomp_flag, "").split(recomp_suffix, 1)[0]).split('_'))
        elif _val.startswith('log_ground_state_'):
            recomp_flag = "log_ground_state_"
            recomp_suffix = ".dat"
            decomposed_list.append(str(_val.replace(recomp_flag, "").split(recomp_suffix, 1)[0]).split('_'))
        elif _val.startswith('state_'):
            recomp_flag = "state_"
            recomp_suffix = ".pkl"
            decomposed_list.append(str(_val.replace(recomp_flag, "").split(recomp_suffix, 1)[0]).split('_'))
        elif _val.startswith('E_psi_M_'):
            recomp_flag = "E_psi_M_"
            recomp_suffix = ".pkl"
            decomposed_list.append(str(_val.replace(recomp_flag, "").split(recomp_suffix, 1)[0]).split('_'))
        else:
            print(_val)

    # add nu entries to the end of each decomposed dat file
    for _i, _val in enumerate(mylist):
        # print(decomposed_list[_i])
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

    sorted_list = sorted(decomposed_list, key=itemgetter(model_index, chi_index, r_index, s_index,
                                                         p_index, q_index, Ly_index))

    # convert the types of entries in the list back into strings
    for _i, _val in enumerate(mylist):
        sorted_list[_i] = [str(x) for x in sorted_list[_i]]

    # remove nu entries and recompose list
    recomposed_list = []
    for _i, _val in enumerate(sorted_list):
        sorted_list[_i] = sorted_list[_i][:len(sorted_list[_i]) - 3]
        recomposed_list.append('_'.join(sorted_list[_i]))

    return recomposed_list


if __name__ == '__main__':

    # get the complete directory list
    complete_dir_list = os.listdir()

    # get the dat/pkl file list
    dat_file_list = []
    suffixes = [".dat", ".pkl"]
    for i, val in enumerate(complete_dir_list):
        if any(x in val for x in suffixes):
            if "Nmax" in val:
                continue
            else:
                dat_file_list.append(val)

    # sort the dat file list
    sorted_dat_file_list = sort_list(dat_file_list)

    with open("list_configurations.out", 'w') as file1:
        for i, val in enumerate(sorted_dat_file_list):
            file1.write(val+"\n")
            print(val)
