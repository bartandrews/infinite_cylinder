import os
import numpy as np

if __name__ == '__main__':

    # initialize lists for r=-4,-3,-2,-1,1,2,3
    r = [[], [], [], [], [], [], [], [], []]

    # populate r=-5 list (5 particles for C=1)
    # r[0].append({"chi": 250, "V": 10, "n": (5, 27), "nphi": (1, 3), "LxMUC": 1, "Ly": 9, "phi_max": 9})

    # populate r=-4 list (8 particles for C=1)
    # r[1].append({"chi": 200, "V": 10, "n": (4, 21), "nphi": (1, 3), "LxMUC": 1, "Ly": 14, "phi_max": 7})

    # populate r=-3 list (6 particles for C=1, 3 particles for C=2)
    # for chi_val in [200, 250]:
    #     r[2].append({"chi": chi_val, "V": 10, "n": (3, 55), "nphi": (3, 5), "LxMUC": 1, "Ly": 11, "phi_max": 11})
    #     r[2].append({"chi": chi_val, "V": 10, "n": (3, 77), "nphi": (4, 7), "LxMUC": 1, "Ly": 11, "phi_max": 11})
    #     r[2].append({"chi": chi_val, "V": 10, "n": (1, 33), "nphi": (5, 9), "LxMUC": 1, "Ly": 11, "phi_max": 11})
    #     r[2].append({"chi": chi_val, "V": 10, "n": (3, 121), "nphi": (6, 11), "LxMUC": 1, "Ly": 11, "phi_max": 11})
    #     r[2].append({"chi": chi_val, "V": 10, "n": (3, 143), "nphi": (7, 13), "LxMUC": 1, "Ly": 11, "phi_max": 11})
    #     r[2].append({"chi": chi_val, "V": 10, "n": (1, 55), "nphi": (8, 15), "LxMUC": 1, "Ly": 11, "phi_max": 11})

    # populate r=-2 list (4 particles for C=1, 4 particles for C=2)
    for chi_val in [100, 150]:
        r[3].append({"chi": chi_val, "V": 10, "n": (1, 44), "nphi": (3, 8), "LxMUC": 1, "Ly": 11, "phi_max": 11})
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 121), "nphi": (4, 11), "LxMUC": 1, "Ly": 11, "phi_max": 11})
        r[3].append({"chi": chi_val, "V": 10, "n": (1, 77), "nphi": (5, 14), "LxMUC": 1, "Ly": 11, "phi_max": 11})
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 187), "nphi": (6, 17), "LxMUC": 1, "Ly": 11, "phi_max": 11})
        r[3].append({"chi": chi_val, "V": 10, "n": (1, 110), "nphi": (7, 20), "LxMUC": 1, "Ly": 11, "phi_max": 11})
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 253), "nphi": (8, 23), "LxMUC": 1, "Ly": 11, "phi_max": 11})

    # populate r=-1 list (2 particles for C=2)
    for chi_val in [100, 150]:
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 40), "nphi": (3, 8), "LxMUC": 1, "Ly": 10, "phi_max": 5})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 55), "nphi": (4, 11), "LxMUC": 1, "Ly": 10, "phi_max": 5})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 70), "nphi": (5, 14), "LxMUC": 1, "Ly": 10, "phi_max": 5})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 85), "nphi": (6, 17), "LxMUC": 1, "Ly": 10, "phi_max": 5})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 100), "nphi": (7, 20), "LxMUC": 1, "Ly": 10, "phi_max": 5})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 15), "nphi": (8, 23), "LxMUC": 1, "Ly": 10, "phi_max": 5})

    # populate r=1 list (2 particles for C=1, 2 particles for C=2)
    for chi_val in [100, 150]:
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 56), "nphi": (3, 8), "LxMUC": 1, "Ly": 14, "phi_max": 7})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 77), "nphi": (4, 11), "LxMUC": 1, "Ly": 14, "phi_max": 7})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 98), "nphi": (5, 14), "LxMUC": 1, "Ly": 14, "phi_max": 7})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 119), "nphi": (6, 17), "LxMUC": 1, "Ly": 14, "phi_max": 7})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 140), "nphi": (7, 20), "LxMUC": 1, "Ly": 14, "phi_max": 7})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 161), "nphi": (8, 23), "LxMUC": 1, "Ly": 14, "phi_max": 7})

    # populate r=2 list (4 particles for C=1, 2 particles for C=2)
    for chi_val in [100, 150]:
        r[6].append({"chi": chi_val, "V": 10, "n": (1, 52), "nphi": (3, 8), "LxMUC": 1, "Ly": 13, "phi_max": 13})
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 143), "nphi": (4, 11), "LxMUC": 1, "Ly": 13, "phi_max": 13})
        r[6].append({"chi": chi_val, "V": 10, "n": (1, 91), "nphi": (5, 14), "LxMUC": 1, "Ly": 13, "phi_max": 13})
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 221), "nphi": (6, 17), "LxMUC": 1, "Ly": 13, "phi_max": 13})
        r[6].append({"chi": chi_val, "V": 10, "n": (1, 130), "nphi": (7, 20), "LxMUC": 1, "Ly": 13, "phi_max": 13})
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 229), "nphi": (8, 23), "LxMUC": 1, "Ly": 13, "phi_max": 13})

    # populate r=3 list (6 particles for C=1)
    # r[7].append({"chi": 250, "V": 10, "n": (1, 7), "nphi": (1, 3), "LxMUC": 1, "Ly": 14, "phi_max": 7})

    # populate r=4 list (4 particles for C=1)
    # r[8].append({"chi": 250, "V": 10, "n": (4, 27), "nphi": (1, 3), "LxMUC": 1, "Ly": 9, "phi_max": 9})

    # pick r index
    ridx = 6

    # define the directory
    if ridx == 0:
        base = "r_-5"
    elif ridx == 1:
        base = "r_-4"
    elif ridx == 2:
        base = "r_-3"
    elif ridx == 3:
        base = "r_-2"
    elif ridx == 4:
        base = "r_-1"
    elif ridx == 5:
        base = "r_1"
    elif ridx == 6:
        base = "r_2"
    elif ridx == 7:
        base = "r_3"
    elif ridx == 8:
        base = "r_4"
    else:
        raise ValueError("ridx has to be between 0 and 8.")

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
