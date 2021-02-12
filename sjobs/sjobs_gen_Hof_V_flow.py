import os
import numpy as np

if __name__ == '__main__':

    # initialize lists for r=-4,-3,-2,-1,1,2,3
    r = [[], [], [], [], [], [], []]

    # populate r=-4 list (nu=4/7, 8 particles)
    r[0].append({"chi": 250, "n": (4, 21), "nphi": (1, 3), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[0].append({"chi": 250, "n": (1, 7), "nphi": (1, 4), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[0].append({"chi": 250, "n": (4, 35), "nphi": (1, 5), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[0].append({"chi": 250, "n": (2, 21), "nphi": (1, 6), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[0].append({"chi": 250, "n": (4, 49), "nphi": (1, 7), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[0].append({"chi": 250, "n": (1, 14), "nphi": (1, 8), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})

    # populate r=-3 list (nu=3/5, 6 particles)
    r[1].append({"chi": 250, "n": (1, 5), "nphi": (1, 3), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[1].append({"chi": 250, "n": (3, 20), "nphi": (1, 4), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[1].append({"chi": 250, "n": (3, 25), "nphi": (1, 5), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[1].append({"chi": 250, "n": (1, 10), "nphi": (1, 6), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[1].append({"chi": 250, "n": (3, 35), "nphi": (1, 7), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[1].append({"chi": 250, "n": (3, 40), "nphi": (1, 8), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})

    # populate r=-2 list (nu=2/3, 4 particles)
    r[2].append({"chi": 250, "n": (2, 9), "nphi": (1, 3), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[2].append({"chi": 250, "n": (1, 6), "nphi": (1, 4), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[2].append({"chi": 250, "n": (2, 15), "nphi": (1, 5), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    # r[2].append({"chi": 250, "n": (1, 9), "nphi": (1, 6), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[2].append({"chi": 250, "n": (2, 21), "nphi": (1, 7), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[2].append({"chi": 250, "n": (1, 12), "nphi": (1, 8), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})

    # populate r=-1 list ()

    # populate r=1 list (nu=1/3, 2 particles)
    r[4].append({"chi": 250, "n": (1, 9), "nphi": (1, 3), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[4].append({"chi": 250, "n": (1, 12), "nphi": (1, 4), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[4].append({"chi": 250, "n": (1, 15), "nphi": (1, 5), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[4].append({"chi": 250, "n": (1, 18), "nphi": (1, 6), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[4].append({"chi": 250, "n": (1, 21), "nphi": (1, 7), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[4].append({"chi": 250, "n": (1, 24), "nphi": (1, 8), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": 1, "V_samp": 41})

    # populate r=2 list (nu=2/5, 4 particles)
    r[5].append({"chi": 250, "n": (2, 15), "nphi": (1, 3), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[5].append({"chi": 250, "n": (1, 10), "nphi": (1, 4), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[5].append({"chi": 250, "n": (2, 25), "nphi": (1, 5), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[5].append({"chi": 250, "n": (1, 15), "nphi": (1, 6), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[5].append({"chi": 250, "n": (2, 35), "nphi": (1, 7), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[5].append({"chi": 250, "n": (1, 20), "nphi": (1, 8), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": 1, "V_samp": 41})

    # populate r=3 list (nu=3/7, 6 particles)
    r[6].append({"chi": 250, "n": (1, 7), "nphi": (1, 3), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[6].append({"chi": 250, "n": (3, 28), "nphi": (1, 4), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[6].append({"chi": 250, "n": (3, 35), "nphi": (1, 5), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[6].append({"chi": 250, "n": (1, 14), "nphi": (1, 6), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[6].append({"chi": 250, "n": (3, 49), "nphi": (1, 7), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})
    r[6].append({"chi": 250, "n": (3, 56), "nphi": (1, 8), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": 1, "V_samp": 41})

    # pick r index
    ridx = 6

    # define the directory
    if ridx == 0:
        base = "r_-4"
    elif ridx == 1:
        base = "r_-3"
    elif ridx == 2:
        base = "r_-2"
    elif ridx == 3:
        base = "r_-1"
    elif ridx == 4:
        base = "r_1"
    elif ridx == 5:
        base = "r_2"
    elif ridx == 6:
        base = "r_3"
    else:
        raise ValueError("ridx has to be between 0 and 6.")

    # loop over configurations
    for j in range(len(r[ridx])):
        file = open(os.path.join(f"Hof_V_flow/{base}", f"chi_{r[ridx][j]['chi']}_"
                                 f"n_{r[ridx][j]['n'][0]}_{r[ridx][j]['n'][1]}_"
                                 f"nphi_{r[ridx][j]['nphi'][0]}_{r[ridx][j]['nphi'][1]}_"
                                 f"LxMUC_{r[ridx][j]['LxMUC']}_Ly_{r[ridx][j]['Ly']}_"
                                 f"Vmin_{r[ridx][j]['V_min']}_Vmax_{r[ridx][j]['V_max']}_"
                                 f"Vsamp_{r[ridx][j]['V_samp']}"), "w+")
        file.write("#!/bin/bash\nexport MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
        file.write(f"srun python code/V_flow.py -path -thr 12 -mod FerHofSqu1 -chi {r[ridx][j]['chi']} "
                   f"-t1 1 -Vtype Coulomb -Vrange 1 "
                   f"-n {r[ridx][j]['n'][0]} {r[ridx][j]['n'][1]} "
                   f"-nphi {r[ridx][j]['nphi'][0]} {r[ridx][j]['nphi'][1]} "
                   f"-LxMUC {r[ridx][j]['LxMUC']} -Ly {r[ridx][j]['Ly']} "
                   f"-V_min {r[ridx][j]['V_min']} -V_max {r[ridx][j]['V_max']} -V_samp {r[ridx][j]['V_samp']}\n")
        file.close()
