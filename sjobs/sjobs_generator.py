import os
import numpy as np

if __name__ == '__main__':

    # define V0
    V0 = 10

    # initialize lists for nu=1/3, 2/5, 3/7
    nu = [[], [], []]

    # populate nu=1/3 list
    nu[0].append({"chi": 1950, "n": (1, 12), "nphi": (1, 4), "Ly": 9})
    nu[0].append({"chi": 2000, "n": (1, 12), "nphi": (1, 4), "Ly": 9})
    nu[0].append({"chi": 950, "n": (1, 15), "nphi": (1, 5), "Ly": 9})
    nu[0].append({"chi": 1000, "n": (1, 15), "nphi": (1, 5), "Ly": 9})
    nu[0].append({"chi": 1950, "n": (1, 15), "nphi": (1, 5), "Ly": 12})
    nu[0].append({"chi": 2000, "n": (1, 15), "nphi": (1, 5), "Ly": 12})
    nu[0].append({"chi": 950, "n": (1, 18), "nphi": (1, 6), "Ly": 9})
    nu[0].append({"chi": 1000, "n": (1, 18), "nphi": (1, 6), "Ly": 9})
    nu[0].append({"chi": 1950, "n": (1, 18), "nphi": (1, 6), "Ly": 12})
    nu[0].append({"chi": 2000, "n": (1, 18), "nphi": (1, 6), "Ly": 12})
    nu[0].append({"chi": 950, "n": (1, 21), "nphi": (1, 7), "Ly": 9})
    nu[0].append({"chi": 1000, "n": (1, 21), "nphi": (1, 7), "Ly": 9})
    nu[0].append({"chi": 1950, "n": (1, 27), "nphi": (2, 9), "Ly": 9})
    nu[0].append({"chi": 2000, "n": (1, 27), "nphi": (2, 9), "Ly": 9})
    nu[0].append({"chi": 950, "n": (1, 33), "nphi": (2, 11), "Ly": 9})
    nu[0].append({"chi": 1000, "n": (1, 33), "nphi": (2, 11), "Ly": 9})

    # populate nu=2/5 list
    nu[1].append({"chi": 1950, "n": (1, 15), "nphi": (1, 6), "Ly": 10})
    nu[1].append({"chi": 2000, "n": (1, 15), "nphi": (1, 6), "Ly": 10})
    nu[1].append({"chi": 1950, "n": (2, 35), "nphi": (1, 7), "Ly": 10})
    nu[1].append({"chi": 2000, "n": (2, 35), "nphi": (1, 7), "Ly": 10})
    nu[1].append({"chi": 1950, "n": (2, 55), "nphi": (2, 11), "Ly": 10})
    nu[1].append({"chi": 2000, "n": (2, 55), "nphi": (2, 11), "Ly": 10})
    nu[1].append({"chi": 1950, "n": (2, 65), "nphi": (2, 13), "Ly": 10})
    nu[1].append({"chi": 2000, "n": (2, 65), "nphi": (2, 13), "Ly": 10})

    # populate nu=3/7 list
    nu[2].append({"chi": 1950, "n": (3, 35), "nphi": (1, 5), "Ly": 14})
    nu[2].append({"chi": 2000, "n": (3, 35), "nphi": (1, 5), "Ly": 14})
    nu[2].append({"chi": 1950, "n": (3, 70), "nphi": (1, 10), "Ly": 14})
    nu[2].append({"chi": 2000, "n": (3, 70), "nphi": (1, 10), "Ly": 14})
    nu[2].append({"chi": 1950, "n": (3, 77), "nphi": (2, 11), "Ly": 14})
    nu[2].append({"chi": 2000, "n": (3, 77), "nphi": (2, 11), "Ly": 14})
    nu[2].append({"chi": 1950, "n": (3, 91), "nphi": (2, 13), "Ly": 14})
    nu[2].append({"chi": 2000, "n": (3, 91), "nphi": (2, 13), "Ly": 14})

    # pick filling factor index (0 for nu=1/3, 1 for nu=2/5, 2 for nu=3/7)
    fidx = 0

    # define the directory
    if fidx == 0:
        base = "nu_1_3"
    elif fidx == 1:
        base = "nu_2_5"
    elif fidx == 2:
        base = "nu_3_7"
    else:
        raise ValueError("fidx has to be between 0 and 2.")

    # loop over interaction range
    for i in np.linspace(1, 3, 11):
        # loop over configurations
        for j in range(len(nu[fidx])):
            V_val = V0 / np.sqrt(2*np.pi*(nu[fidx][j]['nphi'][0]/nu[fidx][j]['nphi'][1]))
            file = open(os.path.join(f"{base}", f"chi_{nu[fidx][j]['chi']}_V_{V_val:.5f}_Vrange_{i:g}_"
                                     f"n_{nu[fidx][j]['n'][0]}_{nu[fidx][j]['n'][1]}_"
                                     f"nphi_{nu[fidx][j]['nphi'][0]}_{nu[fidx][j]['nphi'][1]}_"
                                     f"Ly_{nu[fidx][j]['Ly']}"), "w+")
            file.write("#!/bin/bash\nexport MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
            file.write(f"srun python code/ground_state.py -path -thr 12 -mod FerHofSqu1 -chi {nu[fidx][j]['chi']} "
                       f"-t1 1 -V {V_val:.5f} -Vtype Coulomb -Vrange {i:g} "
                       f"-n {nu[fidx][j]['n'][0]} {nu[fidx][j]['n'][1]} "
                       f"-nphi {nu[fidx][j]['nphi'][0]} {nu[fidx][j]['nphi'][1]} -LxMUC 1 -Ly {nu[fidx][j]['Ly']}\n")
            file.close()
