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
from matplotlib.patches import Polygon

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}')
# matplotlib.verbose.level = 'debug-annoying'

##################
# Ly/lB function #
##################


def LylB_func(nphi, Ly):
    return np.sqrt(2*np.pi*nphi)*Ly


#############################
# line of best fit function #
#############################


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

############################################
# line of best fit function for R^2 > 0.99 #
############################################


def straight_line_of_best_fit(axis, LylB_list, SvN_list, xval=0.25, yval=1.2):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    print("SvN = m*(Ly/lB) + c")
    print(f"(m, m_err, c, c_err) = ({m:.3f}, {m_err:.3f}, {c:.3f}, {c_err:.3f})")
    xvalues = np.linspace(0, max(LylB_list))
    axis.plot(xvalues, m * xvalues + c, '-', c='k', zorder=0)
    axis.text(xval, yval, "$S_\mathrm{{vN}}={gradient:.3f}(L_y/l_\mathrm{{B}})-(\mathbf{{{intercept:.3f}\pm {cerror:.3f}}})$\n$R^2={rsquared:.5f}$".format(
        gradient=m, intercept=abs(c), cerror=c_err, rsquared=r2_value, fontsize=10))

    return m, m_err, c, c_err, r2_value


def line_of_best_fit_values(LylB_list, SvN_list):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    return m, m_err, c, c_err, r2_value

legend_y = 1.75

# define a list of easily-visible markers
markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
           (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
           '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
           '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

if __name__ == '__main__':

    fig = plt.figure(figsize=(6, 6))
    outer_grid = gridspec.GridSpec(1, 1)
    upper_cell = outer_grid[0, 0]
    upper_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, upper_cell, hspace=1.4)

    # nu=1/3 ###########################################################################################################

    ax = plt.subplot(upper_inner_grid[0])

    model="FerHofSqu1"
    filling="nu_1_3"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/short_range_int/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = True
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]

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
        _, _, _, _, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            straight_line_of_best_fit(ax, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), yval=1.3)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    ax.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

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
                if LylB or LylB_outliers:
                    if not error_bars:
                        if LylB:
                            ax.plot(LylB, SvN, '.', marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        if LylB_outliers:
                            ax.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        if LylB_outliers:
                            ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)

    ax.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4, framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax.set_xlim([0, None])
    ax.set_ylim(-1, 3)
    ax.text(10.5, 0, "$R^2>0.99$")
    ax.text(5, -0.6, "$(L_y/l_B)_\mathrm{min}\equiv 8$")
    metal0 = Polygon(((0, -10), (8, -10), (8, 10), (0, 10)), fc=(0, 0, 0, 0.1))
    ax.add_artist(metal0)

    # nu=2/5 ###########################################################################################################

    ax1 = plt.subplot(upper_inner_grid[1])

    model = "FerHofSqu1"
    filling = "nu_2_5"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/out/short_range_int/V_is_10/{model}_Vrange_1_{filling}_accepted.out'
    error_bars = True
    identify_outliers = False
    systematic_points = False
    Ly_min = 0
    LylB_min = 9

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]

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
            if val[3] > LylB_min:
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
                if data_line[3] > LylB_min:
                    LylB_outliers.append(data_line[3])
                    SvN_outliers.append(data_line[4])
                    SvN_outliers_error.append(data_line[5])
            else:
                if data_line[3] > LylB_min:
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
        _, _, _, _, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.99:
            straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), xval=0.20, yval=1.3)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    ax1.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            # print(data_line[3], critical_LylB)
            if LylB_min <= data_line[3] < critical_LylB:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            elif data_line[3] >= critical_LylB:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if LylB or LylB_outliers:
                    if not error_bars:
                        if LylB:
                            ax1.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        if LylB_outliers:
                            ax1.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                markersize=6)
                    else:  # error_bars
                        if LylB:
                            ax1.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                    label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        if LylB_outliers:
                            ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3,
                                    color='k', marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)

    ax1.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax1.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax1.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax1.set_xlim([0, None])
    ax1.set_xticks(np.arange(0, 12, 1))
    ax1.set_ylim(-1, 3)
    ax1.text(9.65, 0, "$R^2>0.99$")
    ax1.text(6.5, -0.5, "$(L_y/l_B)_\mathrm{min}\equiv 9$")
    ax1.scatter(LylB_func(1 / 7, 10), 1.4414, s=100, facecolors='none', edgecolors='r', zorder=5)
    metal0 = Polygon(((0, -10), (9, -10), (9, 10), (0, 10)), fc=(0, 0, 0, 0.1))
    ax1.add_artist(metal0)

    # nu=3/7 ###########################################################################################################

    ax2 = plt.subplot(upper_inner_grid[2])

    model = "FerHofSqu1"
    filling = "nu_3_7"
    file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14.out'
    error_bars = True
    identify_outliers = False
    systematic_points = False
    Ly_min = 0

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]

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
        _, _, _, _, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > 0.943:
            straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), xval=0.29, yval=1.475)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    # ax2.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

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
                        ax2.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        ax2.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                                 markersize=6)
                    else:  # error_bars
                        ax2.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3,
                                     marker=markers[flux_density_index],
                                     label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        ax2.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3,
                                     color='k', marker=markers[flux_density_index], markersize=6)

    ax2.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4,
               framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=9, labelspacing=0.3, columnspacing=0)
    ax2.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax2.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax2.set_xlim([0, None])
    #ax2.set_ylim(-1, 4)
    ax2.text(12.5, 0, "$L_y= 14$")
    ax2.text(6, -0.75, "$(L_y/l_B)_\mathrm{min}\equiv 10$")
    ax2.scatter(LylB_func(1 / 10, 14), 1.9632, s=100, facecolors='none', edgecolors='r', zorder=5)
    metal1 = Polygon(((0, -10), (10, -10), (10, 10), (0, 10)), fc=(0, 0, 0, 0.1))
    ax2.add_artist(metal1)

    # complete plot ####################################################################################################

    # plt.setp(ax.get_xticklabels(), visible=False)
    # ax.tick_params(axis='both', which='major', labelsize=10)

    fig.text(0, 0.93, "(a) $\\nu=1/3$\n error $<0.1\%$", fontsize=12)
    fig.text(0, 0.61, "(b) $\\nu=2/5$\n error $<0.1\%$", fontsize=12)
    fig.text(0, 0.29, "(c) $\\nu=3/7$\n error $\lesssim 3\%$", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TEE/short_range_int_test.png", bbox_inches='tight', dpi=300)
    plt.show()

