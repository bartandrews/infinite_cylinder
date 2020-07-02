# configurations.py --- generate a sorted list of all commands to run in a given Ly/lB interval (execute on dart)
#
# Conditions:
# - checks against the existing pickles and omits these from the list, therefore it is recommended to execute on dart
# - is a pickle of the same configuration already exists, with ANY chi, then it will be omitted from the list
# - the sort is performed for all valid configurations, but we display only those that do not already exist in the pickles directory

import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction as Frac
import math
import os
import fnmatch
from colorama import Fore
from colorama import Style


def cost(p_val, q_val, Lx_val, Ly_val):
    nphi_val = p_val/q_val
    processing_cost = np.exp(nphi_val*20) * q_val * Lx_val * np.exp(Ly_val)
    return processing_cost


def LylB(nphi_val, Ly_val):
    return np.sqrt(2*np.pi*nphi_val)*Ly_val


def min_diff(list_val):  # returns minimum difference between any pair
    n = len(list_val)  # length of the list
    list_val = sorted(list_val)  # sort array in non-decreasing order
    diff = 10 ** 20  # initialize difference as infinite
    for i in range(n - 1):  # find the min diff by comparing adjacent pairs in sorted array
        if list_val[i + 1] - list_val[i] < diff:
            diff = list_val[i + 1] - list_val[i]
    return diff


if __name__ == '__main__':

    model = "FerHofSqu1"  # desired model
    nu = 1 / 3  # desired filling factor
    chi = 1000  # desired max chi (for command labelling)
    Ly_min, Ly_max = 4, 100  # desired domain of Ly such that Ly_min <= Ly <= Ly_max
    LylB_min, LylB_max = 8, 100  # desired range of Ly/lB such that LylB_min < Ly/lB < LylB_max
    LylB_separation = 0.1  # keep all LylB values at least this distance away from each other
    Nmin = 2  # minimum number of particles required in the system

    counter = 0
    for Ly in range(Ly_min, Ly_max+1):
        for p in range(1, 21):
            for q in range(4, 21):
                if math.gcd(p, q) != 1:  # ensure that nphi is a coprime fraction
                    continue
                nphi = p/q
                if (LylB_min/Ly)**2/(2*np.pi) < nphi < np.minimum(0.4, (LylB_max/Ly)**2/(2*np.pi)):
                    for Lx in range(1, 11):
                        if q * Lx * Ly <= 100:  # limit the system size (for memory consumption)
                            if abs(nu*nphi*q*Lx*Ly - int(nu*nphi*q*Lx*Ly)) < 1e-7:  # if number of particles is an integer then accept, otherwise try a larger Lx
                                if int(nu*nphi*q*Lx*Ly) >= Nmin:  # check that there are at least Nmin particles
                                    if counter == 0:
                                        data = np.array([[Lx, Ly, p, q, LylB(nphi, Ly), cost(p, q, Lx, Ly)]])
                                    else:
                                        if all(abs(i - LylB(nphi, Ly)) >= LylB_separation for i in list(data[:, 4])):
                                            data_line = np.array([[Lx, Ly, p, q, LylB(nphi, Ly), cost(p, q, Lx, Ly)]])
                                            data = np.concatenate((data, data_line))
                                        else:
                                            for i in list(data[:, 4]):  # for all the existing LylB values
                                                if abs(i - LylB(nphi, Ly)) < LylB_separation:
                                                    current_config = [Lx, Ly, p, q, LylB(nphi, Ly), cost(p, q, Lx, Ly)]
                                                    existing_config = data[np.where(data[:, 4] == i)].tolist()[0]
                                                    if current_config[5] < existing_config[5]:  # if the cost of the current config is lower than the existing config
                                                        provisional_data = np.copy(data)
                                                        provisional_data[np.where(provisional_data[:, 4] == i)] = np.array([current_config])  # tentatively replace rows
                                                        if abs(min_diff(provisional_data[:, 4])) >= LylB_separation:  # if all points still have the required separation
                                                            data = np.copy(provisional_data)  # actually replace rows
                                                            print("Replaced", existing_config, "with", current_config)
                                                            break
                                    counter += 1
                                    break

    # sort the array by cost
    sorted_array = data[np.argsort(data[:, 5])]

    # normalize the cost column
    sorted_array[:, 5] = sorted_array[:, 5] / sorted_array[0][5]

    # print to the screen the first 50 configurations in ascending cost
    # print to the file the commands for those configurations that do not already exist in pickle_path
    pickle_path = f"/home/bart/PycharmProjects/infinite_cylinder/pickles/ground_state/{model}"
    counter = 0
    nu_frac = Frac(str(nu)).limit_denominator(10)
    with open(f"{model}_nu_{nu_frac.numerator}_{nu_frac.denominator}_allowed.out", "w") as file2:
        print("Lx\tLy\tp\tq\tLylB\tcost")
        with open("configurations.out", 'w') as file:
            for line in sorted_array:
                n = Frac(str(nu*line[2]/line[3])).limit_denominator(100)
                if counter < 20:
                    if any(fnmatch.fnmatch(x, f"*_n_{n.numerator}_{n.denominator}_nphi_{int(line[2])}_{int(line[3])}_LxMUC_{int(line[0])}_Ly_{int(line[1])}*") for x in os.listdir(pickle_path)):
                        print(f"{Fore.GREEN}{int(line[0])}\t{int(line[1])}\t{int(line[2])}\t{int(line[3])}\t{line[4]:.3f}\t{line[5]:.3f}\t{46*(int(line[3])*int(line[1]))/84}{Style.RESET_ALL}")
                    else:
                        print(f"{Fore.RED}{int(line[0])}\t{int(line[1])}\t{int(line[2])}\t{int(line[3])}\t{line[4]:.3f}\t{line[5]:.3f}\t{46*2.25*(int(line[3])*int(line[1]))/84}{Style.RESET_ALL}")
                        for chi_val in [chi-50, chi]:
                            if "Bos" in model:
                                file.write(f"echo python code/ground_state.py -thr 1 -mod {model} -chi {chi_val} -t1 1 -n {n.numerator} {n.denominator} -nphi {int(line[2])} {int(line[3])} -LxMUC {int(line[0])} -Ly {int(line[1])}\n")
                            else:  # "Fer"
                                file.write(f"echo python code/ground_state.py -thr 1 -mod {model} -chi {chi_val} -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n {n.numerator} {n.denominator} -nphi {int(line[2])} {int(line[3])} -LxMUC {int(line[0])} -Ly {int(line[1])}\n")
                    file2.write(f"{int(line[2])}\t{int(line[3])}\t{int(line[1])}\t-\t-\t-\n")
                    counter += 1

    # plot the figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    LylB = [i[4] for i in sorted_array]
    y = range(1, len(LylB)+1)
    ax.scatter(LylB, y, c=y)
    ax.set_xlabel("$L_y/l_B$")
    ax.set_ylabel("$n_\phi$ index (sorted by increasing numerical cost)")
    plt.show()
