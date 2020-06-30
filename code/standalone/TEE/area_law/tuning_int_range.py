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
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
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

if __name__ == '__main__':

    model="FerHofSqu1"
    filling="nu_1_3"
    c_1_3 = []
    c_err_1_3 = []

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_1_{filling}_accepted.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = True

    # set Ly_min
    Ly_min = 0

    fig = plt.figure(figsize=(6, 10))

    outer_grid = gridspec.GridSpec(1, 1)
    upper_cell = outer_grid[0, 0]

    upper_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, upper_cell, hspace=0.7)

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

    ax = plt.subplot(upper_inner_grid[0])

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

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
                            ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)
    ax.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4,
              framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)

    ####################################################################################################################
    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_2_{filling}_accepted.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = True

    # set Ly_min
    Ly_min = 0

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
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
        c_1_3.append(c)
        c_err_1_3.append(c_err)
        if r2value > 0.99:
            straight_line_of_best_fit(ax, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
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
                            ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='g',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        # ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    ####################################################################################################################
    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_3_{filling}_accepted.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = True

    # set Ly_min
    Ly_min = 0

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
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
        if r2value > 0.5:
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
                            ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='r',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        # ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    ax.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax.set_xlim([0, None])
    ax.set_ylim([None, 4])
    ax.axhline(-np.log(np.sqrt(3)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax.text(10, -np.log(np.sqrt(3))-0.1, "Abelian theory", backgroundcolor='white')
    #ax.text(7.3, 3, f"$R_1^2={R1_1_3:.3f}$, $R_2^2={R2_1_3:.3f}$, $R_3^2={R3_1_3:.3f}$")
    ax.text(7.3, 3.3, f"$R_1^2={R1_1_3:.3f}$", c='k')
    ax.text(7.3, 2.65, f"$R_2^2={R2_1_3:.3f}$", c='g')
    ax.text(7.3, 2, f"$R_3^2={R3_1_3:.3f}$", c='r')

    ####################################################################################################################
    ####################################################################################################################

    left, bottom, width, height = [0.21, 0.82, 0.3, 0.05]
    ax3 = fig.add_axes([left, bottom, width, height])

    print(c_1_3, c_err_1_3)

    ax3.errorbar(1, abs(c_1_3[0]), yerr=c_err_1_3[0], ls='none', capsize=3, color='k', marker='.', markersize=6)
    ax3.errorbar(2, abs(c_1_3[1]), yerr=c_err_1_3[1], ls='none', capsize=3, color='g', marker='.', markersize=6)
    ax3.errorbar(3, abs(c_1_3[2]), yerr=c_err_1_3[2], ls='none', capsize=3, color='r', marker='.', markersize=6)
    ax3.axhline(np.log(np.sqrt(3)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax3.set_xlabel("$\kappa$", fontsize=11, labelpad=0)
    ax3.set_ylabel("$\gamma$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"
    c_2_5 = []
    c_err_2_5 = []

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_1_{filling}_accepted.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = False

    # set Ly_min
    Ly_min = 0

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

    ax1 = plt.subplot(upper_inner_grid[1])

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
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
        if r2value > 0.9:
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
            if 9.47 <= data_line[3] < critical_LylB:
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
                            ax1.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='k',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        # ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)
    ax1.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4,
              framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)

    ####################################################################################################################
    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_2_{filling}_total.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = False

    # set Ly_min
    Ly_min = 0

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
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
        c_2_5.append(c)
        c_err_2_5.append(c_err)
        if r2value > 0.5:
            straight_line_of_best_fit(ax1, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
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
                            ax1.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax1.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            ax1.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='g',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        # ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    ####################################################################################################################
    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_2_{filling}_total.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = False

    # set Ly_min
    Ly_min = 0

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
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
        if r2value > 0.5:
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
                            ax1.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax1.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            ax1.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='r',
                                    marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                    markersize=6)
                        # ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    ax1.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax1.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax1.set_xlim([0, None])
    ax1.set_ylim([None, 4])
    ax1.axhline(-np.log(np.sqrt(5)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax1.text(8, -np.log(np.sqrt(5)) - 0.1, "Abelian theory", backgroundcolor='white')
    ax1.text(5.8, 3.3, f"$R_1^2={R1_2_5:.3f}$", c='k')
    ax1.text(5.8, 2.6, f"$R_2^2={R2_2_5:.3f}$", c='g')
    ax1.text(5.8, 1.9, f"$R_3^2={R3_2_5:.3f}$", c='r')

    ####################################################################################################################
    ####################################################################################################################

    left, bottom, width, height = [0.21, 0.5225, 0.3, 0.05]
    ax4 = fig.add_axes([left, bottom, width, height])

    print(c_2_5, c_err_2_5)

    ax4.errorbar(1, abs(c_2_5[0]), yerr=c_err_2_5[0], ls='none', capsize=3, color='k', marker='.', markersize=6)
    ax4.errorbar(2, abs(c_2_5[1]), yerr=c_err_2_5[1], ls='none', capsize=3, color='g', marker='.', markersize=6)
    ax4.errorbar(3, abs(c_2_5[2]), yerr=c_err_2_5[2], ls='none', capsize=3, color='r', marker='.', markersize=6)
    ax4.axhline(np.log(np.sqrt(5)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax4.set_xlabel("$\kappa$", fontsize=11, labelpad=0)
    ax4.set_ylabel("$\gamma$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"
    c_3_7 = []
    c_err_3_7 = []

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14_min.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = False

    # set Ly_min
    Ly_min = 0

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

    ax2 = plt.subplot(upper_inner_grid[2])

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
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
                            ax2.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax2.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            ax2.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='k',
                                     marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                     markersize=6)
                        # ax2.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)
    ax2.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.5, legend_y), handletextpad=0, borderpad=0.4,
               framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)

    ####################################################################################################################
    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14_min.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = False

    # set Ly_min
    Ly_min = 0

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
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
        c_3_7.append(c)
        c_err_3_7.append(c_err)
        if r2value > 0.99:
            straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'g')
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
                            ax2.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax2.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            ax2.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='g',
                                     marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                     markersize=6)
                        # ax2.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    ####################################################################################################################
    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_3_7"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/code/standalone/TEE/area_law/{model}_Vrange_1_{filling}_total_Ly_14_min.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = False

    # set Ly_min
    Ly_min = 0

    ####################################################################################################################

    if not identify_outliers:
        LylB_outlier_values = []
    else:  # identify_outliers
        # BosHofSqu1 nu=1/2
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(2 / 7, 6), LylB_func(4/13, 4), LylB_func(1/4, 6), LylB_func(2/9, 6)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4)]
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6)]
        LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
                               LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
                               LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
                               LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
                               LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
    ####################################################################################################################

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

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # gather the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
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
            straight_line_of_best_fit(ax2, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), 'r')
            c_3_7.append(c)
            c_err_3_7.append(c_err)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    R3_3_7 = r2value

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
                            ax2.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                        # ax2.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                    else:  # error_bars
                        if LylB:
                            ax2.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, color='r',
                                     marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$',
                                     markersize=6)
                        # ax2.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k', marker=markers[flux_density_index], markersize=6)

    ax2.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax2.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax2.set_xlim([0, None])
    ax2.set_ylim([None, 5])
    ax2.axhline(-np.log(np.sqrt(7)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax2.text(12, -np.log(np.sqrt(7)) - 0.1, "Abelian theory", backgroundcolor='white')
    ax2.text(8.5, 4.25, f"$R_1^2={R1_3_7:.3f}$", c='k')
    ax2.text(8.5, 3.5, f"$R_2^2={R2_3_7:.3f}$", c='g')
    ax2.text(8.5, 2.75, f"$R_3^2={R3_3_7:.3f}$", c='r')

    ####################################################################################################################
    ####################################################################################################################

    left, bottom, width, height = [0.21, 0.226, 0.3, 0.05]
    ax5 = fig.add_axes([left, bottom, width, height])

    print(c_3_7, c_err_3_7)

    ax5.errorbar(1, abs(c_3_7[0]), yerr=c_err_3_7[0], ls='none', capsize=3, color='k', marker='.', markersize=6)
    ax5.errorbar(2, abs(c_3_7[1]), yerr=c_err_3_7[1], ls='none', capsize=3, color='g', marker='.', markersize=6)
    ax5.errorbar(3, abs(c_3_7[2]), yerr=c_err_3_7[2], ls='none', capsize=3, color='r', marker='.', markersize=6)
    ax5.axhline(np.log(np.sqrt(7)), color='k', ls='dashed', linewidth=1, zorder=-1, label='theory')
    ax5.set_xlabel("$\kappa$", fontsize=11, labelpad=0)
    ax5.set_ylabel("$\gamma$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################
    ####################################################################################################################

    # plt.setp(ax.get_xticklabels(), visible=False)

    # ax.tick_params(axis='both', which='major', labelsize=10)

    fig.text(0, 0.92, "(a) $\\nu=1/3$\n error $<0.1\%$", fontsize=12)
    fig.text(0, 0.6225, "(b) $\\nu=2/5$\n error $<0.1\%$", fontsize=12)
    fig.text(0, 0.325, "(c) $\\nu=3/7$\n error $<3.5\%$", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TEE/figures/tuning_int_range.png", bbox_inches='tight', dpi=300)
    plt.show()

