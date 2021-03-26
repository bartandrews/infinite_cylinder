import os
import numpy as np

if __name__ == '__main__':

    # initialize lists for r=-5, -4,-3,-2,-1,1,2,3, 4
    r = [[], [], [], [], [], [], [], [], []]

    # populate r=-5 list (nu=5/9)

    # populate r=-4 list (nu=4/7, 8 particles)

    # populate r=-3 list (nu=3/5, 6 particles)
    # for V_val in [10, 1, 0.1]:
    #     r[2].append({"chi": 50, "n": (1, 55), "nphi": (8, 15), "LxMUC": 1, "Ly": 11, "V_min": 0, "V_max": V_val, "V_samp": 41})
    #     r[2].append({"chi": 100, "n": (3, 143), "nphi": (7, 13), "LxMUC": 1, "Ly": 11, "V_min": 0, "V_max": V_val, "V_samp": 41})
    #     r[2].append({"chi": 150, "n": (3, 121), "nphi": (6, 11), "LxMUC": 1, "Ly": 11, "V_min": 0, "V_max": V_val, "V_samp": 41})

    # populate r=-2 list (nu=2/3, 4 particles)
    # for V_val in [10, 1, 0.1]:
    #     r[3].append({"chi": 100, "n": (2, 105), "nphi": (8, 15), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": V_val, "V_samp": 41})
    # for V_val in [10, 1, 0.1]:
    #     r[3].append({"chi": 50, "n": (2, 121), "nphi": (4, 11), "LxMUC": 1, "Ly": 11, "V_min": 0, "V_max": V_val, "V_samp": 41})
    #     r[3].append({"chi": 100, "n": (2, 121), "nphi": (4, 11), "LxMUC": 1, "Ly": 11, "V_min": 0, "V_max": V_val, "V_samp": 41})
    for V_val in [10, 1, 0.1]:
        r[3].append({"chi": 25, "n": (2, 225), "nphi": (4, 15), "LxMUC": 1, "Ly": 15, "V_min": 0, "V_max": V_val, "V_samp": 41})
        r[3].append({"chi": 50, "n": (2, 225), "nphi": (4, 15), "LxMUC": 1, "Ly": 15, "V_min": 0, "V_max": V_val, "V_samp": 41})
        r[3].append({"chi": 25, "n": (2, 361), "nphi": (4, 19), "LxMUC": 1, "Ly": 19, "V_min": 0, "V_max": V_val, "V_samp": 41})

    # populate r=-1 list ()
    # for V_val in [10, 1, 0.1]:
    #     r[4].append({"chi": 25, "n": (1, 27), "nphi": (5, 9), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": V_val, "V_samp": 41})
    #     r[4].append({"chi": 25, "n": (1, 33), "nphi": (6, 11), "LxMUC": 1, "Ly": 6, "V_min": 0, "V_max": V_val, "V_samp": 41})
    # for V_val in [10, 1, 0.1]:
    #     r[4].append({"chi": 50, "n": (1, 85), "nphi": (6, 17), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": V_val, "V_samp": 41})
    #     r[4].append({"chi": 50, "n": (1, 100), "nphi": (7, 20), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": V_val, "V_samp": 41})
    for V_val in [10, 1, 0.1]:
        r[4].append({"chi": 25, "n": (1, 133), "nphi": (5, 19), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": V_val, "V_samp": 41})
        r[4].append({"chi": 25, "n": (1, 161), "nphi": (6, 23), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": V_val, "V_samp": 41})
        r[4].append({"chi": 25, "n": (1, 216), "nphi": (5, 24), "LxMUC": 1, "Ly": 18, "V_min": 0, "V_max": V_val, "V_samp": 41})

    # populate r=1 list (nu=1/3, 2 particles)
    # for V_val in [10, 1, 0.1]:
    #     r[5].append({"chi": 50, "n": (1, 55), "nphi": (6, 11), "LxMUC": 1, "Ly": 10, "V_min": 0, "V_max": V_val, "V_samp": 41})
    # for V_val in [10, 1, 0.1]:
    #     r[5].append({"chi": 50, "n": (1, 98), "nphi": (5, 14), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": V_val, "V_samp": 41})
    #     r[5].append({"chi": 100, "n": (1, 98), "nphi": (5, 14), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": V_val, "V_samp": 41})
    #     r[5].append({"chi": 50, "n": (1, 140), "nphi": (7, 20), "LxMUC": 1, "Ly": 14, "V_min": 0, "V_max": V_val, "V_samp": 41})

    # populate r=2 list (nu=2/5, 4 particles)
    # for V_val in [10, 1, 0.1]:
    #     r[6].append({"chi": 50, "n": (2, 81), "nphi": (5, 9), "LxMUC": 1, "Ly": 9, "V_min": 0, "V_max": V_val, "V_samp": 41})

    # populate r=3 list (nu=3/7, 6 particles)

    # populate r=4 list (nu=4/9)

    # pick r index
    ridx = 4

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
