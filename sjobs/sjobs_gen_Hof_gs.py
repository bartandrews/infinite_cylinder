import os
import numpy as np

if __name__ == '__main__':

    # initialize lists for r=-4,-3,-2,-1,1,2,3
    r = [[], [], [], [], [], [], [], [], []]

    # populate r=-5 list (5 particles for C=1)
    for chi_val in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
        r[0].append({"chi": chi_val, "V": 10, "n": (5, 72), "nphi": (1, 8), "LxMUC": 1, "Ly": 18})

    # populate r=-4 list (8 particles for C=1)
    for chi_val in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
        r[1].append({"chi": chi_val, "V": 10, "n": (2, 21), "nphi": (1, 6), "LxMUC": 1, "Ly": 14})
        r[1].append({"chi": chi_val, "V": 10, "n": (4, 49), "nphi": (1, 7), "LxMUC": 1, "Ly": 14})
        r[1].append({"chi": chi_val, "V": 10, "n": (1, 14), "nphi": (1, 8), "LxMUC": 1, "Ly": 14})

    # populate r=-3 list (6 particles for C=1, 3 particles for C=2)
    for chi_val in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
        r[2].append({"chi": chi_val, "V": 10, "n": (3, 25), "nphi": (1, 5), "LxMUC": 1, "Ly": 10})
        r[2].append({"chi": chi_val, "V": 10, "n": (1, 10), "nphi": (1, 6), "LxMUC": 1, "Ly": 10})
        r[2].append({"chi": chi_val, "V": 10, "n": (3, 35), "nphi": (1, 7), "LxMUC": 1, "Ly": 10})
        r[2].append({"chi": chi_val, "V": 10, "n": (3, 40), "nphi": (1, 8), "LxMUC": 1, "Ly": 10})

    # populate r=-2 list (4 particles for C=1, 4 particles for C=2)
    for chi_val in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 15), "nphi": (1, 5), "LxMUC": 1, "Ly": 6})
        r[3].append({"chi": chi_val, "V": 10, "n": (1, 9), "nphi": (1, 6), "LxMUC": 1, "Ly": 6})
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 21), "nphi": (1, 7), "LxMUC": 1, "Ly": 6})
        r[3].append({"chi": chi_val, "V": 10, "n": (1, 12), "nphi": (1, 8), "LxMUC": 1, "Ly": 6})
    for chi_val in [50, 75, 100, 125, 150, 175, 200, 225, 250]:
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 225), "nphi": (4, 15), "LxMUC": 1, "Ly": 15})
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 361), "nphi": (4, 19), "LxMUC": 1, "Ly": 19})
    for chi_val in [75, 100, 125, 150, 175, 200, 225, 250, 275]:
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 225), "nphi": (4, 15), "LxMUC": 1, "Ly": 15})
    for chi_val in [125, 150, 175, 200, 225, 250, 275, 300, 325]:
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 105), "nphi": (8, 15), "LxMUC": 1, "Ly": 14})
    for chi_val in [175, 200, 225, 250, 275, 300, 325, 350, 375]:
        r[3].append({"chi": chi_val, "V": 10, "n": (2, 105), "nphi": (8, 15), "LxMUC": 1, "Ly": 14})

    # populate r=-1 list (2 particles for C=2)
    for chi_val in [50, 75, 100, 125, 150, 175, 200, 225, 250]:
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 27), "nphi": (5, 9), "LxMUC": 1, "Ly": 6})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 33), "nphi": (6, 11), "LxMUC": 1, "Ly": 6})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 133), "nphi": (5, 19), "LxMUC": 1, "Ly": 14})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 161), "nphi": (6, 23), "LxMUC": 1, "Ly": 14})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 216), "nphi": (5, 24), "LxMUC": 1, "Ly": 18})
    for chi_val in [75, 100, 125, 150, 175, 200, 225, 250, 275]:
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 85), "nphi": (6, 17), "LxMUC": 1, "Ly": 10})
        r[4].append({"chi": chi_val, "V": 10, "n": (1, 100), "nphi": (7, 20), "LxMUC": 1, "Ly": 10})

    # populate r=1 list (2 particles for C=1, 2 particles for C=2)
    for chi_val in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 9), "nphi": (1, 3), "LxMUC": 1, "Ly": 6})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 12), "nphi": (1, 4), "LxMUC": 1, "Ly": 6})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 15), "nphi": (1, 5), "LxMUC": 1, "Ly": 6})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 18), "nphi": (1, 6), "LxMUC": 1, "Ly": 6})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 21), "nphi": (1, 7), "LxMUC": 1, "Ly": 6})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 24), "nphi": (1, 8), "LxMUC": 1, "Ly": 6})
    for chi_val in [75, 100, 125, 150, 175, 200, 225, 250, 275]:
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 55), "nphi": (6, 11), "LxMUC": 1, "Ly": 10})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 98), "nphi": (5, 14), "LxMUC": 1, "Ly": 14})
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 140), "nphi": (7, 20), "LxMUC": 1, "Ly": 14})
    for chi_val in [125, 150, 175, 200, 225, 250, 275, 300, 325]:
        r[5].append({"chi": chi_val, "V": 10, "n": (1, 98), "nphi": (5, 14), "LxMUC": 1, "Ly": 14})

    # populate r=2 list (4 particles for C=1, 2 particles for C=2)
    for chi_val in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 15), "nphi": (1, 3), "LxMUC": 1, "Ly": 10})
        r[6].append({"chi": chi_val, "V": 10, "n": (1, 10), "nphi": (1, 4), "LxMUC": 1, "Ly": 10})
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 25), "nphi": (1, 5), "LxMUC": 1, "Ly": 10})
        r[6].append({"chi": chi_val, "V": 10, "n": (1, 15), "nphi": (1, 6), "LxMUC": 1, "Ly": 10})
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 35), "nphi": (1, 7), "LxMUC": 1, "Ly": 10})
        r[6].append({"chi": chi_val, "V": 10, "n": (1, 20), "nphi": (1, 8), "LxMUC": 1, "Ly": 10})
    for chi_val in [75, 100, 125, 150, 175, 200, 225, 250, 275]:
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 81), "nphi": (5, 9), "LxMUC": 1, "Ly": 9})
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 121), "nphi": (4, 11), "LxMUC": 1, "Ly": 11})
    for chi_val in [125, 150, 175, 200, 225, 250, 275, 300, 325]:
        r[6].append({"chi": chi_val, "V": 10, "n": (2, 121), "nphi": (4, 11), "LxMUC": 1, "Ly": 11})

    # populate r=3 list (6 particles for C=1)
    for chi_val in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
        r[7].append({"chi": chi_val, "V": 10, "n": (3, 35), "nphi": (1, 5), "LxMUC": 1, "Ly": 14})
        r[7].append({"chi": chi_val, "V": 10, "n": (1, 14), "nphi": (1, 6), "LxMUC": 1, "Ly": 14})
        r[7].append({"chi": chi_val, "V": 10, "n": (3, 49), "nphi": (1, 7), "LxMUC": 1, "Ly": 14})
        r[7].append({"chi": chi_val, "V": 10, "n": (3, 56), "nphi": (1, 8), "LxMUC": 1, "Ly": 14})
    for chi_val in [75, 100, 125, 150, 175, 200, 225, 250, 275]:
        r[7].append({"chi": chi_val, "V": 10, "n": (3, 121), "nphi": (6, 11), "LxMUC": 1, "Ly": 11})
    for chi_val in [125, 150, 175, 200, 225, 250, 275, 300, 325]:
        r[7].append({"chi": chi_val, "V": 10, "n": (3, 121), "nphi": (6, 11), "LxMUC": 1, "Ly": 11})
    for chi_val in [175, 200, 225, 250, 275, 300, 325, 350, 375]:
        r[7].append({"chi": chi_val, "V": 10, "n": (3, 121), "nphi": (6, 11), "LxMUC": 1, "Ly": 11})

    # populate r=4 list (4 particles for C=1)
    for chi_val in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
        r[8].append({"chi": chi_val, "V": 10, "n": (1, 18), "nphi": (1, 8), "LxMUC": 1, "Ly": 18})

    # pick r index
    ridx = 8

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
        file = open(os.path.join(f"Hof_gs/{base}", f"chi_{r[ridx][j]['chi']}_V_{r[ridx][j]['V']}_"
                                 f"n_{r[ridx][j]['n'][0]}_{r[ridx][j]['n'][1]}_"
                                 f"nphi_{r[ridx][j]['nphi'][0]}_{r[ridx][j]['nphi'][1]}_"
                                 f"LxMUC_{r[ridx][j]['LxMUC']}_Ly_{r[ridx][j]['Ly']}"), "w+")
        file.write("#!/bin/bash\nexport MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
        file.write(f"srun python code/ground_state.py -path -thr 4 -mod FerHofSqu1 -chi {r[ridx][j]['chi']} "
                   f"-t1 1 -V {r[ridx][j]['V']} -Vtype Coulomb -Vrange 1 "
                   f"-n {r[ridx][j]['n'][0]} {r[ridx][j]['n'][1]} "
                   f"-nphi {r[ridx][j]['nphi'][0]} {r[ridx][j]['nphi'][1]} "
                   f"-LxMUC {r[ridx][j]['LxMUC']} -Ly {r[ridx][j]['Ly']}\n")
        file.close()
