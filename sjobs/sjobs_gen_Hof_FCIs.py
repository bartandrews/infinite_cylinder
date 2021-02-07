import os
import numpy as np

if __name__ == '__main__':

    # initialize lists for r=-3,-2,-1,1,2,3
    r = [[], [], [], [], [], []]

    # populate r=-3 list
    r[0].append({"chi": 300, "V": 10, "n": (3, 25), "nphi": (1, 5), "LxMUC": 1, "Ly": 5, "phi_max": 5})  # 200, +1
    r[0].append({"chi": 100, "V": 10, "n": (1, 33), "nphi": (5, 9), "LxMUC": 1, "Ly": 11, "phi_max": 11})  # +1
    r[0].append({"chi": 100, "V": 1, "n": (3, 77), "nphi": (4, 7), "LxMUC": 1, "Ly": 11, "phi_max": 11})  # basic

    # populate r=-2 list
    r[1].append({"chi": 300, "V": 10, "n": (2, 15), "nphi": (1, 5), "LxMUC": 1, "Ly": 6, "phi_max": 3})  # 200, +1
    r[1].append({"chi": 100, "V": 10, "n": (2, 49), "nphi": (4, 7), "LxMUC": 1, "Ly": 7, "phi_max": 7})  # 50, basic
    r[1].append({"chi": 100, "V": 1, "n": (2, 121), "nphi": (4, 11), "LxMUC": 1, "Ly": 11, "phi_max": 11})  # basic

    # populate r=-1 list
    r[2].append({"chi": 300, "V": 10, "n": (1, 4), "nphi": (1, 4), "LxMUC": 1, "Ly": 6, "phi_max": 1})  # basic
    r[2].append({"chi": 100, "V": 10, "n": (1, 21), "nphi": (4, 7), "LxMUC": 1, "Ly": 9, "phi_max": 3})  # basic
    r[2].append({"chi": 100, "V": 10, "n": (1, 55), "nphi": (4, 11), "LxMUC": 1, "Ly": 10, "phi_max": 5})  # basic
    r[2].append({"chi": 50, "V": 10, "n": (1, 105), "nphi": (4, 15), "LxMUC": 1, "Ly": 14, "phi_max": 7})  # basic

    # populate r=1 list
    r[3].append({"chi": 300, "V": 10, "n": (1, 12), "nphi": (1, 4), "LxMUC": 1, "Ly": 6, "phi_max": 3})  # 50, basic
    r[3].append({"chi": 100, "V": 10, "n": (1, 65), "nphi": (7, 13), "LxMUC": 1, "Ly": 10, "phi_max": 5})  # 50, +3
    r[3].append({"chi": 100, "V": 10, "n": (1, 98), "nphi": (5, 14), "LxMUC": 1, "Ly": 14, "phi_max": 7})  # 50, +1

    # populate r=2 list
    r[4].append({"chi": 300, "V": 10, "n": (1, 10), "nphi": (1, 4), "LxMUC": 1, "Ly": 10, "phi_max": 5})  # 100, basic
    r[4].append({"chi": 100, "V": 1, "n": (2, 63), "nphi": (4, 7), "LxMUC": 1, "Ly": 9, "phi_max": 9})  # basic

    # populate r=3 list
    r[5].append({"chi": 300, "V": 10, "n": (3, 28), "nphi": (1, 4), "LxMUC": 1, "Ly": 7, "phi_max": 7})  # 100, basic
    r[5].append({"chi": 100, "V": 1, "n": (3, 91), "nphi": (4, 7), "LxMUC": 1, "Ly": 13, "phi_max": 13})  # basic

    # pick r index
    ridx = 0

    # define the directory
    if ridx == 0:
        base = "r_-3"
    elif ridx == 1:
        base = "r_-2"
    elif ridx == 2:
        base = "r_-1"
    elif ridx == 3:
        base = "r_1"
    elif ridx == 4:
        base = "r_2"
    elif ridx == 5:
        base = "r_3"
    else:
        raise ValueError("ridx has to be between 0 and 5.")

    # loop over configurations
    for j in range(len(r[ridx])):
        file = open(os.path.join(f"Hof_FCIs/{base}", f"chi_{r[ridx][j]['chi']}_V_{r[ridx][j]['V']}_"
                                 f"n_{r[ridx][j]['n'][0]}_{r[ridx][j]['n'][1]}_"
                                 f"nphi_{r[ridx][j]['nphi'][0]}_{r[ridx][j]['nphi'][1]}_"
                                 f"LxMUC_{r[ridx][j]['LxMUC']}_Ly_{r[ridx][j]['Ly']}_"
                                 f"phi_max_{r[ridx][j]['phi_max']}"), "w+")
        file.write("#!/bin/bash\nexport MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
        file.write(f"srun python code/phi_flow.py -path -thr 12 -mod FerHofSqu1 -chi {r[ridx][j]['chi']} "
                   f"-t1 1 -V {r[ridx][j]['V']} -Vtype Coulomb -Vrange 1 "
                   f"-n {r[ridx][j]['n'][0]} {r[ridx][j]['n'][1]} "
                   f"-nphi {r[ridx][j]['nphi'][0]} {r[ridx][j]['nphi'][1]} "
                   f"-LxMUC {r[ridx][j]['LxMUC']} -Ly {r[ridx][j]['Ly']} "
                   f"-phi_min 0 -phi_max {r[ridx][j]['phi_max']} -phi_samp {r[ridx][j]['phi_max']*10 + 1}\n")
        file.close()
