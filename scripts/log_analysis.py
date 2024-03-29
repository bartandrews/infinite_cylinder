# log_analysis.py --- analyse the log_observables files (execute in directory)
#
# Conditions:
# - tagged files are ignored
# - this script assumes the same number of entries in each file name

# --- python imports
import os
import re
from itertools import groupby
import numpy as np
from colorama import Fore
from colorama import Style
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
        if _val.endswith('.dat'):  # ignore tagged files
            decomposed_list.append(str(_val.replace("log_observables_", "").split(".dat", 1)[0]).split('_'))

    # add nu entries to the end of each decomposed dat file
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
    V_index = decomposed_list[0].index("V") + 1
    Vrange_index = decomposed_list[0].index("V") + 3
    r_index = decomposed_list[0].index("nu") + 1
    s_index = decomposed_list[0].index("nu") + 2
    p_index = decomposed_list[0].index("nphi") + 1
    q_index = decomposed_list[0].index("nphi") + 2
    LxMUC_index = decomposed_list[0].index("LxMUC") + 1
    Ly_index = decomposed_list[0].index("Ly") + 1
    chi_index = decomposed_list[0].index("chi") + 1
    sorted_list = sorted(decomposed_list, key=itemgetter(model_index, r_index, s_index, Vrange_index,
                                                         p_index, q_index, Ly_index, LxMUC_index, V_index, chi_index))

    # convert the types of entries in the list back into strings
    for _i, _val in enumerate(mylist):
        sorted_list[_i] = [str(x) for x in sorted_list[_i]]

    # remove nu entries and recompose list
    recomposed_list = []
    for _i, _val in enumerate(sorted_list):
        sorted_list[_i] = sorted_list[_i][:len(sorted_list[_i]) - 3]
        recomposed_list.append("log_observables_" + '_'.join(sorted_list[_i]) + ".dat")

    return recomposed_list


if __name__ == '__main__':

    # get the complete directory list
    complete_dir_list = os.listdir()

    # get the dat file list
    dat_file_list = []
    for i, val in enumerate(complete_dir_list):
        if ".dat" in val:
            dat_file_list.append(val)

    # sort the dat file list
    sorted_dat_file_list = sort_list(dat_file_list)

    # group by string before the "_chi_*_" (by model)
    model_grouped_list = [list(i) for j, i in groupby(sorted_dat_file_list, key=model_keyf)]

    # subgroup by string after the "_chi_*_" (by system)
    system_grouped_list = []
    for i in range(len(model_grouped_list)):
        system_grouped_list.append(
            [list(i) for j, i in groupby(model_grouped_list[i], key=system_keyf)])

    ####################################################################################################################

    acceptable_data_points, total_data_points, frac_nu_previous, Vrange_previous = 0, 0, 0, 0

    headings = ['model', 'V', 'Vrange', 'nu', 'LxMUC', 'Ly', 'nphi', 'Ly/lB', '2nd_chi', 'max_chi', 'SvN_estimate', 'SvN_error / %', 'status']
    print("{: <11} {: <6} {: <6} {: <6} {: <6} {: <6} {: <6} {: <20} {: <8} {: <8} {: <20} {: <20} {: <6}".format(*headings))

    for i in range(len(system_grouped_list)):  # loop over models
        for j in range(len(system_grouped_list[i])):  # loop over systems
            largest_chi, second_largest_chi = 0, 0

            for k in range(len(system_grouped_list[i][j])):  # search for max chi system
                chi = int(re.split('_chi_(\d+)_', system_grouped_list[i][j][k])[1])
                if chi > largest_chi:
                    max_chi_system = system_grouped_list[i][j][k]
                    largest_chi = chi
            max_chi = int(re.split('_chi_(\d+)_', max_chi_system)[1])

            if len(system_grouped_list[i][j]) > 1:  # if there is more than one chi system
                for k in range(len(system_grouped_list[i][j])):  # search for second max chi system
                    chi = int(re.split('_chi_(\d+)_', system_grouped_list[i][j][k])[1])
                    if second_largest_chi < chi < max_chi:
                        second_max_chi_system = system_grouped_list[i][j][k]
                        second_largest_chi = chi
                second_max_chi = int(re.split('_chi_(\d+)_', second_max_chi_system)[1])
            else:
                second_max_chi = "---"

            # pull system parameters
            debased_dat = str(system_grouped_list[i][j][k].replace("log_observables_", "").split(".dat", 1)[0])
            debased_dat_entries = debased_dat.split('_')
            model = debased_dat_entries[0]
            V = float(debased_dat_entries[debased_dat_entries.index("V") + 1])
            Vrange = debased_dat_entries[debased_dat_entries.index("V") + 3]
            LxMUC = int(debased_dat_entries[debased_dat_entries.index("LxMUC") + 1])
            Ly = int(debased_dat_entries[debased_dat_entries.index("Ly") + 1])
            nn = int(debased_dat_entries[debased_dat_entries.index("n") + 1])
            nd = int(debased_dat_entries[debased_dat_entries.index("n") + 2])
            p = int(debased_dat_entries[debased_dat_entries.index("nphi") + 1])
            q = int(debased_dat_entries[debased_dat_entries.index("nphi") + 2])
            nphi = p / q
            nu = (nn / nd) / nphi
            frac_nu = Frac(str(nu)).limit_denominator(100)
            LylB = np.sqrt(2 * np.pi * nphi) * Ly

            for line in list(open(max_chi_system)):  # find the max SnV
                if "SvN" in line:
                    SvN_line = line.rstrip()
                    max_SvN = float(SvN_line.split()[-1])
                    break

            if len(system_grouped_list[i][j]) > 1:  # if there is more than one chi system
                for line in list(open(second_max_chi_system)):  # find the second max SvN
                    if "SvN" in line:
                        SvN_line = line.rstrip()
                        second_max_SvN = float(SvN_line.split()[-1])
                        break

                lower_SvN = max_SvN
                m = (second_max_SvN - max_SvN) / ((1 / second_max_chi) - (1 / max_chi))
                upper_SvN = max_SvN - m * (1 / max_chi)
                SvN_estimate = (lower_SvN + upper_SvN) / 2
                SvN_error = (upper_SvN - lower_SvN) / 2
                SvN_perc_error = (SvN_error / SvN_estimate) * 100  # percentage error
            else:
                SvN_estimate = max_SvN
                # SvN_error = float('nan')
                SvN_perc_error = float('nan')

            # manually adjust the accept threshold
            if frac_nu.numerator == 2 and frac_nu.denominator == 5:
                accept_threshold = 0.1  # 0.65
            elif frac_nu.numerator == 3 and frac_nu.denominator == 7:
                accept_threshold = 3.5  # 0.1
            else:
                accept_threshold = 0.1

            # compute the status
            if isinstance(SvN_perc_error, float) and abs(SvN_perc_error) < accept_threshold:
                status = f"{Fore.GREEN}OK{Style.RESET_ALL}"
            else:
                status = f"{Fore.RED}ERROR{Style.RESET_ALL}"

            # write to total file, if status is green also write to accepted file
            if Vrange != Vrange_previous or frac_nu != frac_nu_previous:  # if the Vrange or nu is different, open new files
                total_file = open(os.path.join('out', f'{model}_Vrange_{Vrange}_nu_{frac_nu.numerator}_{frac_nu.denominator}_total.out'), 'w')
                accepted_file = open(os.path.join('out', f'{model}_Vrange_{Vrange}_nu_{frac_nu.numerator}_{frac_nu.denominator}_accepted.out'), 'w')
            if Vrange != Vrange_previous:
                Vrange_previous = Vrange
            data_line = f"{p}\t{q}\t{Ly}\t{LylB:.15f}\t{SvN_estimate:.15f}\t{abs(SvN_error):.15f}\t{V:.5f}\n"
            if abs(V-10) < 10e-5:  # if V!=10,  abs(V-10) < 10e-5
                total_file.write(data_line)
                if status == f"{Fore.GREEN}OK{Style.RESET_ALL}":
                    accepted_file.write(data_line)

            # if nu has changed and it's not the last line, print the report and reinitialize stats
            if frac_nu != frac_nu_previous and j != len(system_grouped_list[i]) - 1:
                if frac_nu_previous != 0:
                    # frac_nu_previous = Frac(str(nu_previous)).limit_denominator(100)
                    print(f"Total number of acceptable data points "
                          f"for the nu={frac_nu_previous.numerator:d}/{frac_nu_previous.denominator:d} {model} model"
                          f" = {acceptable_data_points}/{total_data_points}")
                    acceptable_data_points, total_data_points = 0, 0
                frac_nu_previous = frac_nu

            # write to terminal
            data = [model, V, Vrange, nu, LxMUC, Ly, nphi, LylB, second_max_chi, max_chi, SvN_estimate, SvN_perc_error, status]
            print("{: <11} {: <6} {: <6} {: <6} {: <6} {: <6} {: <6} {: <20} {: <8} {: <8} {: <20} {: <20} {: <6}"
                  .format(model, '{:<2.3g}'.format(V), Vrange, '{:d}/{:d}'.format(frac_nu.numerator, frac_nu.denominator), LxMUC, Ly,
                          '{:d}/{:d}'.format(p, q), '{:<10.10g}'.format(LylB), second_max_chi, max_chi,
                          '{:<10.10g}'.format(SvN_estimate), '{:<10.10g}'.format(SvN_perc_error), status))

            # increment the counters
            if status == f"{Fore.GREEN}OK{Style.RESET_ALL}":
                acceptable_data_points += 1
            total_data_points += 1

            # if last line, print the report
            if j == len(system_grouped_list[i]) - 1:
                print(f"Total number of acceptable data points "
                      f"for the nu={frac_nu.numerator:d}/{frac_nu.denominator:d} {model} model"
                      f" = {acceptable_data_points}/{total_data_points}")
