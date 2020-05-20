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


def cost(q_val, Lx_val, Ly_val):
    cost_val = q_val * Lx_val * np.exp(Ly_val)
    return cost_val


def LylB(nphi_val, Ly_val):
    return np.sqrt(2*np.pi*nphi_val)*Ly_val


def construct_pkl_path(_model, _chi, _nn, _nd, _p, _q, _LxMUC, _Ly):

    base_path = f"/home/bart/PycharmProjects/infinite_cylinder/pickles/ground_state/{_model}"
    if model == "BosHofSqu1":
        file_name = f"state_{_model}_chi_*_t1_1_n_{_nn}_{_nd}_nphi_{_p}_{_q}_LxMUC_{_LxMUC}_Ly_{_Ly}.pkl"
    elif model == "FerHofSqu1":
        file_name = f"state_{_model}_chi_*_t1_1_V_10_Coulomb_1_n_{_nn}_{_nd}_nphi_{_p}_{_q}_LxMUC_{_LxMUC}_Ly_{_Ly}.pkl"
    else:
        raise ValueError("Unsupported model for construct_pkl_path function.")

    file_name_old = file_name.replace("state", "E_psi_M")

    return base_path, file_name, file_name_old


if __name__ == '__main__':

    model = "FerHofSqu1"  # desired model (for command labelling)
    nu = 1 / 3  # desired filling factor
    chi = 1000  # desired chi (for command labelling)
    Ly_min, Ly_max = 4, 15  # desired domain of Ly such that Ly_min <= Ly <= Ly_max
    LylB_min, LylB_max = 10, 20  # desired range of Ly/lB such that LylB_min < Ly/lB < LylB_max
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
                        if Lx >= 2 and q*Lx*Ly > 150:  # reject higher Lx values if there is a large associated memory cost
                            break
                        if abs(nu*nphi*q*Lx*Ly - int(nu*nphi*q*Lx*Ly)) < 1e-7:  # if number of particles is an integer then accept, otherwise try a larger Lx
                            if int(nu*nphi*q*Lx*Ly) >= Nmin:  # check that there are at least Nmin particles
                                if counter == 0:
                                    data = np.array([[Lx, Ly, p, q, LylB(nphi, Ly), cost(q, Lx, Ly)]])
                                else:
                                    if all(abs(i - LylB(nphi, Ly)) >= LylB_separation for i in list(data[:, 4])):
                                        data_line = np.array([[Lx, Ly, p, q, LylB(nphi, Ly), cost(q, Lx, Ly)]])
                                        data = np.concatenate((data, data_line))
                                counter += 1
                                break

    # sort the array by cost
    sorted_array = data[np.argsort(data[:, 5])]

    # normalize the cost column
    sorted_array[:, 5] = sorted_array[:, 5] / sorted_array[0][5]

    # print to the screen and file only the configurations that do not already exist
    print("Lx\tLy\tp\tq\tLylB\tcost")
    with open("configurations.out", 'w') as file:
        for line in sorted_array:
            n = Frac(str(nu*line[2]/line[3])).limit_denominator(100)
            pickle_path, pickle_file, pickle_file_old = construct_pkl_path(model, chi, n.numerator, n.denominator, int(line[2]), int(line[3]), int(line[0]), int(line[1]))
            pickle_files = [pickle_file, pickle_file_old]
            if any(fnmatch.fnmatch(x, pickle_file) for x in os.listdir(pickle_path)) or any(fnmatch.fnmatch(y, pickle_file_old) for y in os.listdir(pickle_path)):
                continue
            else:
                print(f"{int(line[0])}\t{int(line[1])}\t{int(line[2])}\t{int(line[3])}\t{line[4]:.3f}\t{line[5]:.3f}")
                for chi_val in [chi-50, chi]:
                    if "Bos" in model:
                        file.write(f"echo python code/ground_state.py -thr 1 -mod {model} -chi {chi_val} -t1 1 -n {n.numerator} {n.denominator} -nphi {int(line[2])} {int(line[3])} -LxMUC {int(line[0])} -Ly {int(line[1])}\n")
                    else:  # "Fer"
                        file.write(f"echo python code/ground_state.py -thr 1 -mod {model} -chi {chi_val} -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n {n.numerator} {n.denominator} -nphi {int(line[2])} {int(line[3])} -LxMUC {int(line[0])} -Ly {int(line[1])}\n")

    # plot the figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    LylB = [i[4] for i in sorted_array]
    y = range(1, len(LylB)+1)
    ax.scatter(LylB, y, c=y)
    ax.set_xlabel("$L_y/l_B$")
    ax.set_ylabel("$n_\phi$ index (sorted by increasing numerical cost)")
    plt.show()
