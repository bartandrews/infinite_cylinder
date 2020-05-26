# list_configs.py --- generate a sorted list of all configurations in a directory (execute in directory)

# --- python imports
import os
from fractions import Fraction as Frac
from operator import itemgetter
import ast


########################################
# tryeval (auto convert types in list) #
########################################


def tryeval(value):
    try:
        value = ast.literal_eval(value)
    except ValueError:
        pass
    return value


#################################
# sort_list (sort the filenames #
#################################


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
            if ".h5" in _val:
                suffix = ".h5"
            elif ".pkl" in _val:
                suffix = ".pkl"
            else:
                raise ValueError("Unknown file extension.")
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
        elif os.path.exists("state_" + '_'.join(sorted_list[_i]) + ".h5"):
            recomposed_list.append("state_" + '_'.join(sorted_list[_i]) + ".h5")
        elif os.path.exists("state_" + '_'.join(sorted_list[_i]) + ".pkl"):  # backward compatibility
            recomposed_list.append("state_" + '_'.join(sorted_list[_i]) + ".pkl")
        elif os.path.exists("E_psi_M_" + '_'.join(sorted_list[_i]) + ".pkl"):  # backward compatibility
            recomposed_list.append("E_psi_M_" + '_'.join(sorted_list[_i]) + ".pkl")
        else:
            raise ValueError("Missing file corresponding to the configuration: ", '_'.join(sorted_list[_i]))

    return recomposed_list


if __name__ == '__main__':

    # get the complete directory list
    complete_dir_list = os.listdir()

    # get the (dat, h5, pkl) file list
    file_list = []
    suffixes = [".dat", ".h5", ".pkl"]
    for i, val in enumerate(complete_dir_list):
        if any(x in val for x in suffixes):
            if "Nmax" in val:  # do not consider MR states
                continue
            else:
                file_list.append(val)

    # sort the file list
    sorted_file_list = sort_list(file_list)

    with open("list_configurations.out", 'w') as file1:
        for i, val in enumerate(sorted_file_list):
            file1.write(val+"\n")
            print(val)
