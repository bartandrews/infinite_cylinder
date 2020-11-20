import os
import numpy as np

if __name__ == '__main__':

    # define V0
    V0 = 10

    # initialize lists for nu=1/3, 2/5, 3/7
    nu = [[], [], []]

    # populate nu=2/5 list
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 35), "nphi": (2, 7), "Ly": 5})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (2, 35), "nphi": (2, 7), "Ly": 5})
    nu[1].append({"chi": 400, "Vrange": 1, "n": (2, 35), "nphi": (2, 7), "Ly": 10})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (2, 35), "nphi": (2, 7), "Ly": 10})
    nu[1].append({"chi": 200, "Vrange": 1, "n": (2, 35), "nphi": (2, 7), "Ly": 15})
    nu[1].append({"chi": 850, "Vrange": 1, "n": (2, 35), "nphi": (2, 7), "Ly": 15})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 45), "nphi": (2, 9), "Ly": 5})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (2, 45), "nphi": (2, 9), "Ly": 5})
    nu[1].append({"chi": 1950, "Vrange": 1, "n": (2, 45), "nphi": (2, 9), "Ly": 10})
    nu[1].append({"chi": 2000, "Vrange": 1, "n": (2, 45), "nphi": (2, 9), "Ly": 10})
    nu[1].append({"chi": 200, "Vrange": 1, "n": (2, 45), "nphi": (2, 9), "Ly": 15})
    nu[1].append({"chi": 1950, "Vrange": 1, "n": (2, 55), "nphi": (2, 11), "Ly": 10})
    nu[1].append({"chi": 2000, "Vrange": 1, "n": (2, 55), "nphi": (2, 11), "Ly": 10})
    nu[1].append({"chi": 850, "Vrange": 1, "n": (2, 65), "nphi": (2, 13), "Ly": 10})
    nu[1].append({"chi": 900, "Vrange": 1, "n": (2, 65), "nphi": (2, 13), "Ly": 10})
    nu[1].append({"chi": 1950, "Vrange": 1, "n": (1, 20), "nphi": (3, 8), "Ly": 10})
    nu[1].append({"chi": 2000, "Vrange": 1, "n": (1, 20), "nphi": (3, 8), "Ly": 10})
    nu[1].append({"chi": 1950, "Vrange": 1, "n": (1, 25), "nphi": (3, 10), "Ly": 10})
    nu[1].append({"chi": 2000, "Vrange": 1, "n": (1, 25), "nphi": (3, 10), "Ly": 10})
    nu[1].append({"chi": 850, "Vrange": 1, "n": (1, 55), "nphi": (3, 11), "Ly": 10})
    nu[1].append({"chi": 900, "Vrange": 1, "n": (1, 55), "nphi": (3, 11), "Ly": 10})
    nu[1].append({"chi": 1950, "Vrange": 1, "n": (2, 65), "nphi": (5, 13), "Ly": 6})
    nu[1].append({"chi": 2000, "Vrange": 1, "n": (2, 65), "nphi": (5, 13), "Ly": 6})
    nu[1].append({"chi": 1550, "Vrange": 1, "n": (2, 65), "nphi": (5, 13), "Ly": 7})
    nu[1].append({"chi": 1600, "Vrange": 1, "n": (2, 65), "nphi": (5, 13), "Ly": 7})
    nu[1].append({"chi": 400, "Vrange": 1, "n": (1, 35), "nphi": (5, 14), "Ly": 5})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (1, 35), "nphi": (5, 14), "Ly": 5})
    nu[1].append({"chi": 1950, "Vrange": 1, "n": (1, 35), "nphi": (5, 14), "Ly": 6})
    nu[1].append({"chi": 2000, "Vrange": 1, "n": (1, 35), "nphi": (5, 14), "Ly": 6})
    nu[1].append({"chi": 1950, "Vrange": 1, "n": (1, 35), "nphi": (5, 14), "Ly": 7})
    nu[1].append({"chi": 2000, "Vrange": 1, "n": (1, 35), "nphi": (5, 14), "Ly": 7})
    nu[1].append({"chi": 200, "Vrange": 1, "n": (1, 35), "nphi": (5, 14), "Ly": 10})
    nu[1].append({"chi": 400, "Vrange": 1, "n": (1, 35), "nphi": (5, 14), "Ly": 10})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 40), "nphi": (5, 16), "Ly": 5})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 40), "nphi": (5, 16), "Ly": 6})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 40), "nphi": (5, 16), "Ly": 7})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (1, 40), "nphi": (5, 16), "Ly": 8})
    nu[1].append({"chi": 950, "Vrange": 1, "n": (1, 40), "nphi": (5, 16), "Ly": 8})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 40), "nphi": (5, 16), "Ly": 9})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (1, 40), "nphi": (5, 16), "Ly": 9})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 40), "nphi": (5, 16), "Ly": 10})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 85), "nphi": (5, 17), "Ly": 6})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 85), "nphi": (5, 17), "Ly": 7})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 85), "nphi": (5, 17), "Ly": 8})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 85), "nphi": (5, 17), "Ly": 10})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (2, 85), "nphi": (5, 17), "Ly": 10})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 45), "nphi": (5, 18), "Ly": 6})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 45), "nphi": (5, 18), "Ly": 7})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 45), "nphi": (5, 18), "Ly": 8})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 45), "nphi": (5, 18), "Ly": 9})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (1, 45), "nphi": (5, 18), "Ly": 9})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (1, 45), "nphi": (5, 18), "Ly": 10})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (1, 45), "nphi": (5, 18), "Ly": 10})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 95), "nphi": (5, 19), "Ly": 5})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 95), "nphi": (5, 19), "Ly": 6})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 95), "nphi": (5, 19), "Ly": 7})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 95), "nphi": (5, 19), "Ly": 8})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 95), "nphi": (5, 19), "Ly": 10})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (2, 95), "nphi": (5, 19), "Ly": 10})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 105), "nphi": (5, 21), "Ly": 5})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (2, 105), "nphi": (5, 21), "Ly": 5})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 105), "nphi": (5, 21), "Ly": 6})
    nu[1].append({"chi": 800, "Vrange": 1, "n": (2, 105), "nphi": (5, 21), "Ly": 6})
    nu[1].append({"chi": 750, "Vrange": 1, "n": (2, 105), "nphi": (5, 21), "Ly": 7})
    nu[1].append({"chi": 1950, "Vrange": 2, "n": (2, 55), "nphi": (2, 11), "Ly": 10})
    nu[1].append({"chi": 2000, "Vrange": 2, "n": (2, 55), "nphi": (2, 11), "Ly": 10})
    nu[1].append({"chi": 2950, "Vrange": 3, "n": (2, 55), "nphi": (2, 11), "Ly": 10})
    nu[1].append({"chi": 3000, "Vrange": 3, "n": (2, 55), "nphi": (2, 11), "Ly": 10})
    nu[1].append({"chi": 1950, "Vrange": 2, "n": (2, 65), "nphi": (2, 13), "Ly": 10})
    nu[1].append({"chi": 2000, "Vrange": 2, "n": (2, 65), "nphi": (2, 13), "Ly": 10})
    nu[1].append({"chi": 1950, "Vrange": 3, "n": (2, 65), "nphi": (2, 13), "Ly": 10})
    nu[1].append({"chi": 2000, "Vrange": 3, "n": (2, 65), "nphi": (2, 13), "Ly": 10})

    # populate nu=3/7 list
    nu[2].append({"chi": 1950, "Vrange": 1, "n": (1, 21), "nphi": (2, 9), "Ly": 7})
    nu[2].append({"chi": 2000, "Vrange": 1, "n": (1, 21), "nphi": (2, 9), "Ly": 7})
    nu[2].append({"chi": 1950, "Vrange": 1, "n": (3, 77), "nphi": (2, 11), "Ly": 14})
    nu[2].append({"chi": 2000, "Vrange": 1, "n": (3, 77), "nphi": (2, 11), "Ly": 14})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 91), "nphi": (2, 13), "Ly": 7})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 91), "nphi": (2, 13), "Ly": 7})
    nu[2].append({"chi": 1950, "Vrange": 1, "n": (3, 91), "nphi": (2, 13), "Ly": 14})
    nu[2].append({"chi": 2000, "Vrange": 1, "n": (3, 91), "nphi": (2, 13), "Ly": 14})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (1, 35), "nphi": (2, 15), "Ly": 7})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (1, 35), "nphi": (2, 15), "Ly": 7})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 133), "nphi": (2, 19), "Ly": 7})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 133), "nphi": (2, 19), "Ly": 7})
    nu[2].append({"chi": 1950, "Vrange": 1, "n": (3, 56), "nphi": (3, 8), "Ly": 7})
    nu[2].append({"chi": 2000, "Vrange": 1, "n": (3, 56), "nphi": (3, 8), "Ly": 7})
    nu[2].append({"chi": 3250, "Vrange": 1, "n": (3, 70), "nphi": (3, 10), "Ly": 7})
    nu[2].append({"chi": 3300, "Vrange": 1, "n": (3, 70), "nphi": (3, 10), "Ly": 7})
    nu[2].append({"chi": 2400, "Vrange": 1, "n": (3, 77), "nphi": (3, 11), "Ly": 7})
    nu[2].append({"chi": 2450, "Vrange": 1, "n": (3, 77), "nphi": (3, 11), "Ly": 7})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 91), "nphi": (3, 13), "Ly": 7})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 91), "nphi": (3, 13), "Ly": 7})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 98), "nphi": (3, 14), "Ly": 7})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 98), "nphi": (3, 14), "Ly": 7})
    nu[2].append({"chi": 1950, "Vrange": 1, "n": (3, 77), "nphi": (4, 11), "Ly": 7})
    nu[2].append({"chi": 2000, "Vrange": 1, "n": (3, 77), "nphi": (4, 11), "Ly": 7})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 91), "nphi": (4, 13), "Ly": 7})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (1, 42), "nphi": (5, 18), "Ly": 7})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (1, 42), "nphi": (5, 18), "Ly": 7})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 133), "nphi": (5, 19), "Ly": 7})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 133), "nphi": (5, 19), "Ly": 7})
    nu[2].append({"chi": 1950, "Vrange": 1, "n": (3, 119), "nphi": (6, 17), "Ly": 7})
    nu[2].append({"chi": 2000, "Vrange": 1, "n": (3, 119), "nphi": (6, 17), "Ly": 7})
    nu[2].append({"chi": 1950, "Vrange": 1, "n": (3, 133), "nphi": (6, 19), "Ly": 7})
    nu[2].append({"chi": 2000, "Vrange": 1, "n": (3, 133), "nphi": (6, 19), "Ly": 7})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 4})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 4})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 5})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 5})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 6})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 6})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 8})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 8})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 9})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 9})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (1, 42), "nphi": (7, 18), "Ly": 10})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 133), "nphi": (7, 19), "Ly": 4})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 133), "nphi": (7, 19), "Ly": 4})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 133), "nphi": (7, 19), "Ly": 5})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 133), "nphi": (7, 19), "Ly": 5})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 133), "nphi": (7, 19), "Ly": 6})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 133), "nphi": (7, 19), "Ly": 6})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 140), "nphi": (7, 20), "Ly": 4})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 140), "nphi": (7, 20), "Ly": 4})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 140), "nphi": (7, 20), "Ly": 6})
    nu[2].append({"chi": 800, "Vrange": 1, "n": (3, 140), "nphi": (7, 20), "Ly": 6})
    nu[2].append({"chi": 750, "Vrange": 1, "n": (3, 140), "nphi": (7, 20), "Ly": 8})
    nu[2].append({"chi": 1950, "Vrange": 2, "n": (3, 77), "nphi": (2, 11), "Ly": 14})
    nu[2].append({"chi": 2000, "Vrange": 2, "n": (3, 77), "nphi": (2, 11), "Ly": 14})
    nu[2].append({"chi": 1950, "Vrange": 3, "n": (3, 77), "nphi": (2, 11), "Ly": 14})
    nu[2].append({"chi": 2000, "Vrange": 3, "n": (3, 77), "nphi": (2, 11), "Ly": 14})
    nu[2].append({"chi": 1950, "Vrange": 2, "n": (3, 91), "nphi": (2, 13), "Ly": 14})
    nu[2].append({"chi": 2000, "Vrange": 2, "n": (3, 91), "nphi": (2, 13), "Ly": 14})
    nu[2].append({"chi": 1950, "Vrange": 3, "n": (3, 91), "nphi": (2, 13), "Ly": 14})
    nu[2].append({"chi": 2000, "Vrange": 3, "n": (3, 91), "nphi": (2, 13), "Ly": 14})

    # pick filling factor index (0 for nu=1/3, 1 for nu=2/5, 2 for nu=3/7)
    fidx = 2

    # define the directory
    if fidx == 0:
        base = "nu_1_3"
    elif fidx == 1:
        base = "nu_2_5"
    elif fidx == 2:
        base = "nu_3_7"
    else:
        raise ValueError("fidx has to be between 0 and 2.")

    # loop over configurations
    for j in range(len(nu[fidx])):
        V_val = V0
        file = open(os.path.join(f"fixed_n/{base}", f"chi_{nu[fidx][j]['chi']}_V_{V_val:g}_Vrange_{nu[fidx][j]['Vrange']:g}_"
                                 f"n_{nu[fidx][j]['n'][0]}_{nu[fidx][j]['n'][1]}_"
                                 f"nphi_{nu[fidx][j]['nphi'][0]}_{nu[fidx][j]['nphi'][1]}_"
                                 f"Ly_{nu[fidx][j]['Ly']}"), "w+")
        file.write("#!/bin/bash\nexport MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
        file.write(f"srun python code/ground_state.py -path -thr 12 -mod FerHofSqu1 -chi {nu[fidx][j]['chi']} "
                   f"-t1 1 -V {V_val:g} -Vtype Coulomb -Vrange {nu[fidx][j]['Vrange']:g} "
                   f"-n {nu[fidx][j]['n'][0]} {nu[fidx][j]['n'][1]} "
                   f"-nphi {nu[fidx][j]['nphi'][0]} {nu[fidx][j]['nphi'][1]} -LxMUC 1 -Ly {nu[fidx][j]['Ly']}\n")
        file.close()
