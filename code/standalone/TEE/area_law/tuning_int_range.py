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

    # print("SvN = m*(Ly/lB) + c")
    # print(f"(m, m_err, c, c_err) = ({m:.3f}, {m_err:.3f}, {c:.3f}, {c_err:.3f})")
    # xvalues = np.arange(max(LylB_list) + 1)
    # ax.plot(xvalues, m * xvalues + c, '-', c='k', zorder=0)
    # fig.text(0.4, 0.2, "$S_\mathrm{{vN}}={gradient:.3f}(L_y/l_\mathrm{{B}})+({intercept:.3f}\pm {cerror:.3f})$\n$R^2={rsquared:.10f}$".format(
    #     gradient=m, intercept=c, cerror=c_err, rsquared=r2_value))

    return m, m_err, c, c_err, r2_value


def straight_line_of_best_fit(axis, LylB_list, SvN_list, color='k'):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    print("SvN = m*(Ly/lB) + c")
    print(f"(m, m_err, c, c_err) = ({m:.3f}, {m_err:.3f}, {c:.3f}, {c_err:.3f})")
    xvalues = np.linspace(0, max(LylB_list))
    axis.plot(xvalues, m * xvalues + c, '-', c=color, zorder=0)
    # axis.text(0.25, 1.8, "$S_\mathrm{{vN}}={gradient:.3f}(L_y/l_\mathrm{{B}})+({intercept:.3f}\pm {cerror:.3f})$\n$R^2={rsquared:.5f}$".format(
    #     gradient=m, intercept=c, cerror=c_err, rsquared=r2_value, fontsize=10))

    return m, m_err, c, c_err, r2_value


def line_of_best_fit_values(LylB_list, SvN_list):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    return m, m_err, c, c_err, r2_value


legend_y = 1.35

# define a list of easily-visible markers
markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
           (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
           '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
           '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

if __name__ == '__main__':

    fig = plt.figure(figsize=(6, 10))
    outer_grid = gridspec.GridSpec(1, 1)
    upper_cell = outer_grid[0, 0]
    upper_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, upper_cell, hspace=0.7)

    # nu=1/3, Vrange=1 #################################################################################################

    ax = plt.subplot(upper_inner_grid[0])
    c_1_3 = []
    c_err_1_3 = []

    model="FerHofSqu1"
    filling="nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
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
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0]/a[1])]

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
            straight_line_of_best_fit(ax, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R1_1_3 = r2value

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if 8.1 <= data_line[3] < critical_LylB:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB != [] or LylB_outliers != []:
                    if not error_bars:
                        if LylB:
                            ax.plot(LylB, SvN, '.', marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            flux_density_index += 1
                            if flux_density_index == 5 or flux_density_index == 6:
                                flux_density_index += 2
                            print("flux_density_index_1 = ", flux_density_index)
                            ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)
    ax.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4,
              framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)

    # nu=1/3, Vrange=1.2 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.2_{filling}_accepted.out'
    error_bars = True
    identify_outliers = True
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.98:
            # straight_line_of_best_fit(ax, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R1p2_1_3 = r2value

    # nu=1/3, Vrange=1.4 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.4_{filling}_accepted.out'
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

    R1p4_1_3 = r2value

    # nu=1/3, Vrange=1.6 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.6_{filling}_accepted.out'
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

    R1p6_1_3 = r2value

    # nu=1/3, Vrange=1.8 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.8_{filling}_accepted.out'
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

    R1p8_1_3 = r2value

    # nu=1/3, Vrange=2 #################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_10/{model}_Vrange_2_{filling}_accepted.out'
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
            straight_line_of_best_fit(ax, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R2_1_3 = r2value

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if 8 <= data_line[3] < critical_LylB:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB != [] or LylB_outliers != []:
                    if not error_bars:
                        if LylB:
                            ax.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            flux_density_index += 1
                            if flux_density_index == 5 or flux_density_index == 6:
                                flux_density_index += 2
                            print("flux_density_index_2 = ", flux_density_index)
                            ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='g',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        # ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    # nu=1/3, Vrange=2.2 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2.2_{filling}_accepted.out'
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

    R2p2_1_3 = r2value

    # nu=1/3, Vrange=2.4 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2.4_{filling}_accepted.out'
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

    R2p4_1_3 = r2value

    # nu=1/3, Vrange=2.6 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2.6_{filling}_accepted.out'
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

    R2p6_1_3 = r2value

    # nu=1/3, Vrange=2.8 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2.8_{filling}_accepted.out'
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

    R2p8_1_3 = r2value

    # nu=1/3, Vrange=3 #################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_3_{filling}_total.out'
    error_bars = True
    identify_outliers = True
    systematic_points = True
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 4, 9), LylB_func(1 / 5, 12)]

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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.97:
            straight_line_of_best_fit(ax, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'r')
            c_1_3.append(c)
            c_err_1_3.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R3_1_3 = r2value

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or 8 <= data_line[3] < critical_LylB:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB != [] or LylB_outliers != []:
                    if not error_bars:
                        if LylB:
                            ax.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            # five points
                            # if flux_density_index == 0:  # 4th point
                            #     flux_density_index = 2
                            # elif flux_density_index == 1:  # 2nd point
                            #     flux_density_index = 3
                            # elif flux_density_index == 2:  # 1st point
                            #     flux_density_index = 4
                            # elif flux_density_index == 3:  # 5th point
                            #     flux_density_index = 7
                            # elif flux_density_index == 4:  # 3rd point
                            #     flux_density_index = 8
                            # eight points
                            flux_density_index += 1
                            if flux_density_index == 5 or flux_density_index == 6:
                                flux_density_index += 2
                            print("flux_density_index_3 = ", flux_density_index)
                            ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='r',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        if LylB_outliers:
                            if 11 < LylB_outliers[0] < 12:
                                flux_density_index = 1
                            elif 13 < LylB_outliers[0] < 14:
                                flux_density_index = 2
                            print("flux_density_index_outlier = ", flux_density_index)
                            ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, markerfacecolor='w', markeredgecolor='r', ecolor='r', marker=markers[flux_density_index], markersize=6)

    ax.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax.set_xlim([0, None])
    ax.set_ylim([-1, 4])
    ax.axhline(-np.log(np.sqrt(3)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax.text(10, -1.2*np.log(np.sqrt(3)), "Abelian theory", backgroundcolor='white')
    #ax.text(7.3, 3, f"$R_1^2={R1_1_3:.3f}$, $R_2^2={R2_1_3:.3f}$, $R_3^2={R3_1_3:.3f}$")
    ax.text(7.5, 3.4, f"$R_1^2={R1_1_3:.3f}$", c='k')
    ax.text(7.5, 2.75, f"$R_2^2={R2_1_3:.3f}$", c='g')
    ax.text(7.5, 2.1, f"$R_3^2={R3_1_3:.3f}$", c='r')

    # nu=1/3 subplot ###################################################################################################

    left, bottom, width, height = [0.22, 0.82, 0.3, 0.05]
    ax3 = fig.add_axes([left, bottom, width, height])

    print(c_1_3, c_err_1_3)

    ax3.errorbar(1, abs(c_1_3[0]), yerr=c_err_1_3[0], ls='none', capsize=3, color='k', marker='.', markersize=6)
    ax3.errorbar(1.2, abs(c_1_3[1]), yerr=c_err_1_3[1], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax3.errorbar(1.4, abs(c_1_3[2]), yerr=c_err_1_3[2], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax3.errorbar(1.6, abs(c_1_3[3]), yerr=c_err_1_3[3], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax3.errorbar(1.8, abs(c_1_3[4]), yerr=c_err_1_3[4], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax3.errorbar(2, abs(c_1_3[5]), yerr=c_err_1_3[5], ls='none', capsize=3, color='g', marker='.', markersize=6)
    ax3.errorbar(2.2, abs(c_1_3[6]), yerr=c_err_1_3[6], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax3.errorbar(2.4, abs(c_1_3[7]), yerr=c_err_1_3[7], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax3.errorbar(2.6, abs(c_1_3[8]), yerr=c_err_1_3[8], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax3.errorbar(2.8, abs(c_1_3[9]), yerr=c_err_1_3[9], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax3.errorbar(3, abs(c_1_3[10]), yerr=c_err_1_3[10], ls='none', capsize=3, color='r', marker='.', markersize=6)
    ax3.axhline(np.log(np.sqrt(3)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax3.set_xlabel("$\kappa$", fontsize=11, labelpad=0)
    ax3.set_ylabel("$\gamma$", fontsize=11)

    # nu=2/5, Vrange=1 #################################################################################################

    ax1 = plt.subplot(upper_inner_grid[1])
    c_2_5 = []
    c_err_2_5 = []

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1_{filling}_accepted.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.99:
            straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R1_2_5 = r2value
    print(sorted_plot_data)

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or 9.47 <= data_line[3] < critical_LylB:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB != [] or LylB_outliers != []:
                    if not error_bars:
                        if LylB:
                            ax1.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax1.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            if flux_density_index == 0:
                                flux_density_index = 2
                            elif flux_density_index == 1:
                                flux_density_index = 3
                            elif flux_density_index == 2:
                                flux_density_index = 13
                            elif flux_density_index == 3:
                                flux_density_index = 14
                            print("flux_density_index_4 = ", flux_density_index)
                            ax1.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='k',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        # ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)
    ax1.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4,
              framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)

    # nu=2/5, Vrange=1.2 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.2_{filling}_accepted.out'
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

    R1p2_2_5 = r2value

    # nu=2/5, Vrange=1.4 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.4_{filling}_accepted.out'
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

    R1p4_2_5 = r2value

    # nu=2/5, Vrange=1.6 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.6_{filling}_accepted.out'
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

    R1p6_2_5 = r2value

    # nu=2/5, Vrange=1.8 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.8_{filling}_accepted.out'
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

    R1p8_2_5 = r2value

    # nu=2/5, Vrange=2 #################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2_{filling}_accepted.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.93:
            straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R2_2_5 = r2value

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or 8 <= data_line[3] < critical_LylB:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB != [] or LylB_outliers != []:
                    if not error_bars:
                        if LylB:
                            ax1.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax1.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            if flux_density_index == 0:
                                flux_density_index = 2
                            elif flux_density_index == 1:
                                flux_density_index = 3
                            elif flux_density_index == 2:
                                flux_density_index = 13
                            elif flux_density_index == 3:
                                flux_density_index = 14
                            print("flux_density_index_5 = ", flux_density_index)
                            ax1.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='g',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        # ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    # nu=2/5, Vrange=2.2 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2.2_{filling}_accepted.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.9:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R2p2_2_5 = r2value
    print("2.2 : ", c, c_err)

    # nu=2/5, Vrange=2.4 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2.4_{filling}_accepted.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.9:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R2p4_2_5 = r2value
    print("2.4 : ", c, c_err)

    # nu=2/5, Vrange=2.6 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2.6_{filling}_accepted.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.9:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R2p6_2_5 = r2value
    print("2.6 : ", c, c_err)

    # nu=2/5, Vrange=2.8 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_2.8_{filling}_accepted.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.9:
            # straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R2p8_2_5 = r2value
    print("2.8 : ", c, c_err)

    # nu=2/5, Vrange=3 #################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_3_{filling}_total.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.9:
            straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'r')
            c_2_5.append(c)
            c_err_2_5.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R3_2_5 = r2value

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or 8 <= data_line[3] < critical_LylB:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB != [] or LylB_outliers != []:
                    if not error_bars:
                        if LylB:
                            ax1.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax1.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            if flux_density_index == 0:
                                flux_density_index = 2
                            elif flux_density_index == 1:
                                flux_density_index = 3
                            elif flux_density_index == 2:
                                flux_density_index = 13
                            elif flux_density_index == 3:
                                flux_density_index = 14
                            print("flux_density_index_6 = ", flux_density_index)
                            ax1.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='r',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        if LylB_outliers:
                            flux_density_index = 13
                            print("flux_density_index_outlier = ", flux_density_index)
                            ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3,
                                         markerfacecolor='w', markeredgecolor='r', marker=markers[flux_density_index],
                                         markersize=6)
                        # ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    ax1.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax1.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax1.set_xlim([0, None])
    ax1.set_ylim([None, 4])
    ax1.axhline(-np.log(np.sqrt(5)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax1.text(8, -1.2*np.log(np.sqrt(5)), "Abelian theory", backgroundcolor='white')
    ax1.text(5.5, 3.2, f"$R_1^2={R1_2_5:.3f}$", c='k')
    ax1.text(5.5, 2.4, f"$R_2^2={R2_2_5:.3f}$", c='g')
    ax1.text(5.5, 1.6, f"$R_3^2={R3_2_5:.3f}$", c='r')

    # nu=2/5 subplot ###################################################################################################

    left, bottom, width, height = [0.19, 0.5225, 0.3, 0.05]
    ax4 = fig.add_axes([left, bottom, width, height])

    print(c_2_5, c_err_2_5)

    ax4.errorbar(1, abs(c_2_5[0]), yerr=c_err_2_5[0], ls='none', capsize=3, color='k', marker='.', markersize=6)
    ax4.errorbar(1.2, abs(c_2_5[1]), yerr=c_err_2_5[1], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax4.errorbar(1.4, abs(c_2_5[2]), yerr=c_err_2_5[2], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax4.errorbar(1.6, abs(c_2_5[3]), yerr=c_err_2_5[3], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax4.errorbar(1.8, abs(c_2_5[4]), yerr=c_err_2_5[4], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax4.errorbar(2, abs(c_2_5[5]), yerr=c_err_2_5[5], ls='none', capsize=3, color='g', marker='.', markersize=6)
    ax4.errorbar(2.2, abs(c_2_5[6]), yerr=c_err_2_5[6], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax4.errorbar(2.4, abs(c_2_5[7]), yerr=c_err_2_5[7], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax4.errorbar(2.6, abs(c_2_5[8]), yerr=c_err_2_5[8], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax4.errorbar(2.8, abs(c_2_5[9]), yerr=c_err_2_5[9], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax4.errorbar(3, abs(c_2_5[10]), yerr=c_err_2_5[10], ls='none', capsize=3, color='r', marker='.', markersize=6)
    ax4.axhline(np.log(np.sqrt(5)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax4.set_xlabel("$\kappa$", fontsize=11, labelpad=0)
    ax4.set_ylabel("$\gamma$", fontsize=11)

    # nu=3/7, Vrange=1 #################################################################################################

    ax2 = plt.subplot(upper_inner_grid[2])
    c_3_7 = []
    c_err_3_7 = []

    model = "FerHofSqu1"
    filling = "nu_3_7"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_10/{model}_Vrange_1_{filling}_total.out'
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
            # print(data_line)
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                # print(data_line)
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > 10:
                    print(data_line)
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.99:
            straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R1_3_7 = r2value

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or 8 <= data_line[3] < critical_LylB or data_line[2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB and data_line[3] > 10:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB != [] or LylB_outliers != []:
                    if not error_bars:
                        if LylB:
                            ax2.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax2.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            if flux_density_index == 1:  # 1st point
                                flux_density_index = 0
                            elif flux_density_index == 5:  # 2nd point
                                flux_density_index = 4
                            elif flux_density_index == 12:  # 3rd point
                                flux_density_index = 6
                            elif flux_density_index == 13:  # 4th point
                                flux_density_index = 7
                            print("flux_density_index_7 = ", flux_density_index)
                            ax2.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='k',
                                     marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                     markersize=6)
                        # ax2.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)
    ax2.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4,
               framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)

    # nu=3/7, Vrange=1.2 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.2_{filling}_total.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.99:
            # straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R1p2_3_7 = r2value

    # nu=3/7, Vrange=1.4 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.4_{filling}_total.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.98:
            # straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    print("len(c_3_7) = ", len(c_3_7))
    R1p4_3_7 = r2value

    # nu=3/7, Vrange=1.6 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.6_{filling}_total.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.99:
            # straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R1p6_3_7 = r2value

    # nu=3/7, Vrange=1.8 ###############################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_not_10/{model}_Vrange_1.8_{filling}_total.out'
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
        print(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print(r2value)
        if r2value > 0.99:
            # straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R1p8_3_7 = r2value

    # nu=3/7, Vrange=2 #################################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/tuning_int_range/V_is_10/{model}_Vrange_2_{filling}_total.out'
    error_bars = True
    identify_outliers = False
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
        print("this = ", sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        print("this = ", r2value)
        if r2value > 0.98:
            straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R2_3_7 = r2value

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or 8 <= data_line[3] < critical_LylB:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB != [] or LylB_outliers != []:
                    if not error_bars:
                        if LylB:
                            ax2.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax2.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            if flux_density_index == 0:  # 1st point
                                flux_density_index = 0
                            elif flux_density_index == 1:  # 2nd point
                                flux_density_index = 4
                            elif flux_density_index == 2:  # 3rd point
                                flux_density_index = 6
                            elif flux_density_index == 3:  # 4th point
                                flux_density_index = 7
                            print("flux_density_index_8 = ", flux_density_index)
                            ax2.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='g',
                                     marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                     markersize=6)
                        # ax2.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    # nu=3/7, Vrange=3 #################################################################################################

    # model = "FerHofSqu1"
    # filling = "nu_3_7"
    # file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_3_{filling}_total.out'
    # error_bars = True
    # identify_outliers = False
    # systematic_points = False
    # Ly_min = 14
    #
    # if not identify_outliers:
    #     LylB_outlier_values = []
    # else:  # identify_outliers
    #     LylB_outlier_values = [LylB_func(1 / 6, 14), LylB_func(1 / 7, 14), LylB_func(1 / 8, 14), LylB_func(1 / 11, 14)]
    #
    # # append data from file to a list
    # data = []
    # with open(file, 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     line_count = 0
    #     for row in plots:
    #         data.append(int(row[0]))
    #         data.append(int(row[1]))
    #         data.append(int(row[2]))
    #         data.append(float(row[3]))
    #         data.append(float(row[4]))
    #         data.append(float(row[5]))
    #         line_count += 1
    #
    # # group the list into lines of the file
    # i = 0
    # grouped_data = []
    # while i < len(data):
    #     grouped_data.append(data[i:i + 6])
    #     i += 6
    #
    # # delete lines of data such that only systematic points are left for Ly/lB > 8
    # if systematic_points:
    #     file2 = f'/home/bart/PycharmProjects/infinite_cylinder/scripts/{model}_{filling}_allowed.out'
    #     sys_data = []
    #     with open(file2, 'r') as csvfile:
    #         sys_plots = csv.reader(csvfile, delimiter='\t')
    #         line_count = 0
    #         for row in sys_plots:
    #             sys_data.append(int(row[0]))
    #             sys_data.append(int(row[1]))
    #             sys_data.append(int(row[2]))
    #             line_count += 1
    #     # group the list into lines of the file
    #     i = 0
    #     sys_grouped_data = []
    #     while i < len(sys_data):
    #         sys_grouped_data.append(sys_data[i:i + 3])
    #         i += 3
    #     del_indices = []
    #     for i, val in enumerate(grouped_data):
    #         if val[3] > 8:
    #             if any(val[:3] == allowed for allowed in sys_grouped_data):
    #                 continue
    #             else:
    #                 del_indices.append(i)
    #
    #     new_grouped_data = [i for j, i in enumerate(grouped_data) if j not in del_indices]
    #     grouped_data = new_grouped_data
    #
    # # group data by flux density
    # flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0] / a[1])]
    #
    # # gather the data by flux density
    # for flux_density_index in range(len(flux_grouped_data)):
    #     LylB, SvN, SvN_error = [], [], []
    #     LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
    #     for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
    #         if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
    #             if data_line[3] > 8:
    #                 LylB_outliers.append(data_line[3])
    #                 SvN_outliers.append(data_line[4])
    #                 SvN_outliers_error.append(data_line[5])
    #         else:
    #             if data_line[3] > 8:
    #                 LylB.append(data_line[3])
    #                 SvN.append(data_line[4])
    #                 SvN_error.append(data_line[5])
    #
    # # plot the line of best fit
    # LylB = []
    # SvN = []
    # for i, val in enumerate(grouped_data):
    #     if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
    #         continue
    #     else:
    #         LylB.append(val[3])
    #         SvN.append(val[4])
    #
    # plot_data = np.column_stack((LylB, SvN))
    # sorted_plot_data = plot_data[plot_data[:, 0].argsort()]
    # length = len(sorted_plot_data)
    #
    # for i in range(length):
    #     _, _, c, c_err, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
    #     if r2value > 0.8:
    #         straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'r')
    #         c_3_7.append(c)
    #         c_err_3_7.append(c_err)
    #         critical_LylB = sorted_plot_data[0][0]
    #         break
    #     sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)
    #
    # R3_3_7 = r2value
    #
    # # plot the data by flux density
    # for flux_density_index in range(len(flux_grouped_data)):
    #     LylB, SvN, SvN_error = [], [], []
    #     LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
    #     for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
    #         # print(data_line[3], critical_LylB)
    #         if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or 8 <= data_line[3] < critical_LylB:
    #             LylB_outliers.append(data_line[3])
    #             SvN_outliers.append(data_line[4])
    #             SvN_outliers_error.append(data_line[5])
    #         elif data_line[3] >= critical_LylB:
    #             LylB.append(data_line[3])
    #             SvN.append(data_line[4])
    #             SvN_error.append(data_line[5])
    #         if i == len(flux_grouped_data[flux_density_index]) - 1:
    #             if LylB != [] or LylB_outliers != []:
    #                 if not error_bars:
    #                     if LylB:
    #                         ax2.plot(LylB, SvN, '.', marker=markers[flux_density_index],
    #                              label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
    #                     # ax2.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
    #                 else:  # error_bars
    #                     if LylB:
    #                         if flux_density_index == 0:  # 1st point
    #                             flux_density_index = 0
    #                         elif flux_density_index == 1:  # 2nd point
    #                             flux_density_index = 4
    #                         elif flux_density_index == 2:  # 3rd point
    #                             flux_density_index = 6
    #                         elif flux_density_index == 3:  # 4th point
    #                             flux_density_index = 7
    #                         print("flux_density_index_9 = ", flux_density_index)
    #                         ax2.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='r',
    #                                  marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
    #                                  markersize=6)
    #                     # ax2.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)
    #
    ax2.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax2.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax2.set_xlim([0, None])
    ax2.set_ylim([None, 5])
    ax2.axhline(-np.log(np.sqrt(7)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax2.text(12, -1.2*np.log(np.sqrt(7)), "Abelian theory", backgroundcolor='white')
    ax2.text(8.1, 4.1, f"$R_1^2={R1_3_7:.3f}$", c='k')
    ax2.text(8.1, 3.25, f"$R_2^2={R2_3_7:.3f}$", c='g')
    # #ax2.text(8.1, 2.4, f"$R_3^2={R3_3_7:.3f}$", c='r')

    # nu=3/7 subplot ###################################################################################################

    left, bottom, width, height = [0.19, 0.226, 0.3, 0.05]
    ax5 = fig.add_axes([left, bottom, width, height])

    print(c_3_7, c_err_3_7)

    ax5.errorbar(1, abs(c_3_7[0]), yerr=c_err_3_7[0], ls='none', capsize=3, color='k', marker='.', markersize=6)
    ax5.errorbar(1.2, abs(c_3_7[1]), yerr=c_err_3_7[1], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax5.errorbar(1.4, abs(c_3_7[2]), yerr=c_err_3_7[2], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax5.errorbar(1.6, abs(c_3_7[3]), yerr=c_err_3_7[3], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax5.errorbar(1.8, abs(c_3_7[4]), yerr=c_err_3_7[4], ls='none', capsize=3, color='b', marker='.', markersize=6)
    ax5.errorbar(2, abs(c_3_7[5]), yerr=c_err_3_7[5], ls='none', capsize=3, color='g', marker='.', markersize=6)
    # ax5.errorbar(3, abs(c_3_7[2]), yerr=c_err_3_7[2], ls='none', capsize=3, color='r', marker='.', markersize=6)
    ax5.axhline(np.log(np.sqrt(7)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax5.set_xlabel("$\kappa$", fontsize=11, labelpad=-8)
    ax5.set_xticks([1, 2])
    ax5.set_ylabel("$\gamma$", fontsize=11)

    # complete plot ####################################################################################################

    # plt.setp(ax.get_xticklabels(), visible=False)
    # ax.tick_params(axis='both', which='major', labelsize=10)

    fig.text(0, 0.92, "(a) $\\nu=1/3$\n error $<0.1\%$", fontsize=12)
    fig.text(0, 0.6225, "(b) $\\nu=2/5$\n error $<0.1\%$", fontsize=12)
    fig.text(0, 0.325, "(c) $\\nu=3/7$\n error $\lesssim 3\%$", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TEE/tuning_int_range_new.png", bbox_inches='tight', dpi=300)
    plt.show()

