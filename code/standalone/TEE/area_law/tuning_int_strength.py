import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import groupby
# import sys
import math
from scipy import stats
from uncertainties import unumpy, ufloat
import matplotlib.gridspec as gridspec
import os
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}')
# matplotlib.verbose.level = 'debug-annoying'


def LylB_func(nphi, Ly):
    return np.sqrt(2*np.pi*nphi)*Ly


def line_of_best_fit(LylB_list, SvN_list):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    return m, m_err, c, c_err, r2_value


def diagnostic_plot(title, x, y):

    # scatter plot
    plt.scatter(x, y)

    # line of best fit
    parameters, cov = np.polyfit(x, y, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(x, y)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value * r_value
    xvalues = np.linspace(0, max(x))
    plt.plot(np.array(xvalues), m * np.array(xvalues) + c, '-', c='k', zorder=0)
    plt.text(0.25, 0.9*max(y),
              "$S_\mathrm{{vN}}={gradient:.3f}(L_y/l_\mathrm{{B}})+({intercept:.3f}\pm {cerror:.3f})$\n$R^2={rsquared:.5f}$"
             .format(gradient=m, intercept=c, cerror=c_err, rsquared=r2_value, fontsize=10))

    # formatting
    plt.title(title, fontsize=12)
    plt.xlabel("$L_y/l_B$", fontsize=11)
    plt.xlim(0)
    plt.ylabel("$S_\mathrm{vN}$", fontsize=11)
    plt.show()


# define a list of easily-visible markers
markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
           (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
           '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
           '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']


if __name__ == '__main__':

    fig = plt.figure(figsize=(6, 4))
    outer_grid = gridspec.GridSpec(1, 1)
    upper_cell = outer_grid[0, 0]
    upper_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, upper_cell, hspace=0)

    diagnos_plot = "nu_3_7_V_50"

    # nu=1/3, V=10 #####################################################################################################

    c_1_3 = []
    c_err_1_3 = []

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = True
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 9), LylB_func(1 / 5, 12)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_1_3_V_10":
        diagnostic_plot("$\\nu=1/3$, $V_0=10$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R10_1_3 = r2value

    # nu=1/3, V=20 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = True
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 9), LylB_func(1 / 5, 12)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_1_3_V_20":
        diagnostic_plot("$\\nu=1/3$, $V_0=20$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R20_1_3 = r2value

    # nu=1/3, V=30 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = True
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_1_3_V_30":
        diagnostic_plot("$\\nu=1/3$, $V_0=30$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R30_1_3 = r2value

    # nu=1/3, V=40 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = True
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_1_3_V_40":
        diagnostic_plot("$\\nu=1/3$, $V_0=40$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R40_1_3 = r2value

    # nu=1/3, V=50 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = True
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            # straight_line_of_best_fit(ax, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_1_3_V_50":
        diagnostic_plot("$\\nu=1/3$, $V_0=50$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R50_1_3 = r2value

    # nu=2/5, V=10 #####################################################################################################

    c_2_5 = []
    c_err_2_5 = []

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(2 / 11, 10)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_2_5_V_10":
        diagnostic_plot("$\\nu=2/5$, $V_0=10$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R1p2_2_5 = r2value

    # nu=2/5, V=20 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(2 / 11, 10)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_2_5_V_20":
        diagnostic_plot("$\\nu=2/5$, $V_0=20$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R20_2_5 = r2value

    # nu=2/5, V=30 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_2_5_V_30":
        diagnostic_plot("$\\nu=2/5$, $V_0=30$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R30_2_5 = r2value

    # nu=2/5, V=40 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_2_5_V_40":
        diagnostic_plot("$\\nu=2/5$, $V_0=40$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R40_2_5 = r2value

    # nu=2/5, V=50 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                if data_line[3] > 8:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 8:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_2_5_V_50":
        diagnostic_plot("$\\nu=2/5$, $V_0=50$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R50_2_5 = r2value

    # nu=3/7, V=10 #####################################################################################################

    c_3_7 = []
    c_err_3_7 = []

    model = "FerHofSqu1"
    filling = "nu_3_7"
    # file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14.out'
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = True
    systematic_points = False
    Ly_min = 14

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 6, 14), LylB_func(1 / 7, 14), LylB_func(1 / 8, 14), LylB_func(1 / 11, 14)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 10:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_3_7_V_10":
        diagnostic_plot("$\\nu=3/7$, $V_0=10$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R10_3_7 = r2value

    # nu=3/7, V=20 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    # file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14.out'
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = True
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 6, 14), LylB_func(1 / 7, 14), LylB_func(1 / 8, 14), LylB_func(1 / 11, 14)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 10:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_3_7_V_20":
        diagnostic_plot("$\\nu=3/7$, $V_0=20$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R20_3_7 = r2value

    # nu=3/7, V=30 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    # file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14.out'
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = True
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 6, 14), LylB_func(1 / 7, 14), LylB_func(1 / 8, 14), LylB_func(1 / 11, 14)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 10:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_3_7_V_30":
        diagnostic_plot("$\\nu=3/7$, $V_0=30$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R30_3_7 = r2value

    # nu=3/7, V=40 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    # file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14.out'
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = True
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 6, 14), LylB_func(1 / 7, 14), LylB_func(1 / 8, 14), LylB_func(1 / 11, 14)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 10:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:  # 0.943
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_3_7_V_40":
        diagnostic_plot("$\\nu=3/7$, $V_0=40$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R40_3_7 = r2value

    # nu=3/7, V=50 #####################################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    # file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14.out'
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_strength/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = True
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 6, 14), LylB_func(1 / 7, 14), LylB_func(1 / 8, 14), LylB_func(1 / 11, 14)]

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(int(row[2]))
            data.append(float(row[3]))
            data.append(float(row[4]))
            data.append(float(row[5]))
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 6])
        i += 6

    # delete lines of data such that only systematic points are left for Ly/lB > 8
    if systematic_points:
        file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
        sys_data = []
        with open(file2, 'r') as csvfile:
            sys_plots = csv.reader(csvfile, delimiter='\t')
            line_count = 0
            for row in sys_plots:
                sys_data.append(int(row[0]))
                sys_data.append(int(row[1]))
                sys_data.append(int(row[2]))
                line_count += 1
        # group the list into lines of the file
        i = 0
        sys_grouped_data = []
        while i < len(sys_data):
            sys_grouped_data.append(sys_data[i:i + 3])
            i += 3
        del_indices = []
        for i, val in enumerate(grouped_data):
            if val[3] > 8:
                if any(val[:3] == allowed for allowed in sys_grouped_data):
                    continue
                else:
                    del_indices.append(i)

        new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
        grouped_data = new_grouped_data

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 10:
                    LylB.append(data_line[3])
                    SvN.append(data_line[4])
                    SvN_error.append(data_line[5])

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    plot_data = np.column_stack((LylB, SvN))
    sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    length = len(sorted_plot_data)

    for i in range(length):
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:  # 0.943
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    if diagnos_plot == "nu_3_7_V_50":
        diagnostic_plot("$\\nu=3/7$, $V_0=50$", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())

    R50_3_7 = r2value

    # complete plot ####################################################################################################

    ax = plt.subplot(upper_inner_grid[0])
    ax1 = plt.subplot(upper_inner_grid[1], sharex=ax)
    ax2 = plt.subplot(upper_inner_grid[2], sharex=ax1)

    ax.tick_params('x', direction='in', bottom=True)
    ax1.tick_params('x', direction='in', bottom=True)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)

    V = [10, 20, 30, 40, 50]
    # print("c_1_3 = ", c_1_3)
    # print("c_err_1_3 = ", c_err_1_3)
    # print("c_2_5 = ", c_2_5)
    # print("c_err_2_5 = ", c_err_2_5)
    # print("c_3_7 = ", c_3_7)
    # print("c_err_3_7 = ", c_err_3_7)

    gamma_1_3 = [-x for x in c_1_3]
    gamma_2_5 = [-x for x in c_2_5]
    gamma_3_7 = [-x for x in c_3_7]
    ax.errorbar(V, gamma_1_3, yerr=c_err_1_3, label="$\\nu=1/3$", ls='none', capsize=3, color='C0', marker=markers[0], markersize=6)
    ax1.errorbar(V, gamma_2_5, yerr=c_err_2_5, label="$\\nu=2/5$", ls='none', capsize=3, color='C1', marker=markers[1], markersize=6)
    ax2.errorbar(V, gamma_3_7, yerr=c_err_3_7, label="$\\nu=3/7$", ls='none', capsize=3, color='C2', marker=markers[2], markersize=6)
    ax2.set_xlabel("$V_0$", fontsize=11)
    ax.set_ylabel("$\gamma$", fontsize=11)
    ax1.set_ylabel("$\gamma$", fontsize=11)
    ax2.set_ylabel("$\gamma$", fontsize=11)
    ax.legend(loc='upper right', handletextpad=0, borderpad=0.4,
              framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax1.legend(loc='upper right', handletextpad=0, borderpad=0.4,
              framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax2.legend(loc='upper right', handletextpad=0, borderpad=0.4,
              framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)

    ax.axhline(np.log(np.sqrt(3)), color='k', ls='dashed', linewidth=1, zorder=2)
    ax1.axhline(np.log(np.sqrt(5)), color='k', ls='dashed', linewidth=1, zorder=2)
    ax2.axhline(np.log(np.sqrt(7)), color='k', ls='dashed', linewidth=1, zorder=2)

    ax2.set_xticks(V)

    # plt.savefig("/home/bart/Documents/papers/TEE/tuning_int_strength.png", bbox_inches='tight', dpi=300)
    # plt.show()
