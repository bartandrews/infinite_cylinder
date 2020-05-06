import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction as Frac
import math


def cost(q_val, Lx_val, Ly_val):
    cost_val = q_val * Lx_val * np.exp(Ly_val)
    return cost_val


def LylB(nphi_val, Ly_val):
    return np.sqrt(2*np.pi*nphi_val)*Ly_val


if __name__ == '__main__':

    model = "BosHofSqu1"  # desired model (for command labelling)
    nu = 1 / 2  # desired filling factor
    chi = 800  # desired chi (for command labelling)
    Ly_min, Ly_max = 3, 15  # desired domain of Ly such that Ly_min <= Ly <= Ly_max
    LylB_min, LylB_max = 10, 15  # desired range of Ly/lB such that LylB_min < Ly/lB < LylB_max
    LylB_separation = 0.5  # keep all LylB values at least this distance away from each other
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

    # print to the screen and file
    print("Lx\tLy\tp\tq\tLylB\tcost")
    with open("configurations.out", 'w') as file:
        for line in sorted_array:
            print(f"{int(line[0])}\t{int(line[1])}\t{int(line[2])}\t{int(line[3])}\t{line[4]:.3f}\t{line[5]:.3f}")
            n = Frac(str(nu*line[2]/line[3])).limit_denominator(100)
            for chi_val in [chi - 50, chi]:
                if "Bos" in model:
                    file.write(f"echo python code/ground_state.py -thr 1 -mod {model} -chi {chi_val} "
                               f"-t1 1 -n {n.numerator} {n.denominator} -nphi {int(line[2])} {int(line[3])} "
                               f"-LxMUC {int(line[0])} -Ly {int(line[1])}\n")
                else:  # "Fer"
                    file.write(f"echo python code/ground_state.py -thr 1 -mod {model} -chi {chi_val} "
                               f"-t1 1 -V 10 -Vtype \"Coulomb\" -Vrange 1 -n {n.numerator} {n.denominator} "
                               f"-nphi {int(line[2])} {int(line[3])} -LxMUC {int(line[0])} -Ly {int(line[1])}\n")

    # plot the figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    LylB = [i[4] for i in sorted_array]
    y = range(1, len(LylB)+1)
    ax.scatter(LylB, y, c=y)
    ax.set_xlabel("$L_y/l_B$")
    ax.set_ylabel("$n_\phi$ index (sorted by increasing numerical cost)")
    plt.show()
