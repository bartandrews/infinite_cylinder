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


def straight_line_of_best_fit(axis, LylB_list, SvN_list, xpos=0.25, ypos=1.5):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    print("SvN = m*(Ly/lB) + c")
    print(f"(m, m_err, c, c_err) = ({m:.3f}, {m_err:.3f}, {c:.3f}, {c_err:.3f})")
    xvalues = np.linspace(0, max(LylB_list))
    axis.plot(xvalues, m * xvalues + c, '-', c='k', zorder=0)
    axis.text(xpos, ypos, "$S_\mathrm{{vN}}={gradient:.3f}(L_y/l_\mathrm{{B}})+({intercept:.3f}\pm {cerror:.3f})$\n$R^2={rsquared:.5f}$".format(
        gradient=m, intercept=c, cerror=c_err, rsquared=r2_value, fontsize=10))

    return m, m_err, c, c_err, r2_value


def line_of_best_fit_values(LylB_list, SvN_list):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    return m, m_err, c, c_err, r2_value


if __name__ == '__main__':

    model="BosHofSqu1"
    filling="nu_1_2"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_{filling}_accepted.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

    # plot only the systematically collected points for Ly/lB > 8?
    systematic_points = False

    # set Ly_min
    Ly_min = 0

    fig = plt.figure(figsize=(13.75, 15))

    outer_grid = gridspec.GridSpec(3, 1, hspace=0.5)
    upper_cell = outer_grid[0, 0]
    middle_cell = outer_grid[1, 0]
    lower_cell = outer_grid[2, 0]

    upper_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 2, upper_cell, hspace=0)
    middle_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 2, middle_cell, hspace=0)
    lower_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 2, lower_cell, hspace=0)

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
    ax.tick_params('x', direction='in', bottom=True)

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

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if not error_bars:
                    ax.plot(LylB, SvN, '.', marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                else:  # error_bars
                    ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index], label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r', marker=markers[flux_density_index], markersize=6)

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
            straight_line_of_best_fit(ax, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    ax.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

    ax.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.53, 2), handletextpad=0, borderpad=0.4, framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax.set_xlim([0, None])
    ax.set_ylim(-1, 3)

    ####################################################################################################################

    ax1 = plt.subplot(upper_inner_grid[2], sharex=ax)
    ax1.tick_params('x', direction='in', bottom=True)

    # plot the data by flux density
    LylB_total = []
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB_total.append(data_line[3])
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])

    clist = []
    r2list = []

    xvalues = sorted(LylB_total, key=float)[:len(LylB_total)-2]

    # print(xvalues)

    for LylB_min in xvalues:
        # plot the line of best fit
        LylB = []
        SvN = []
        for i, val in enumerate(grouped_data):
            if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
                continue
            else:
                if val[3] >= LylB_min:
                    LylB.append(val[3])
                    SvN.append(val[4])
        m, m_err, c, c_err, r2_value = line_of_best_fit_values(LylB, SvN)
        clist.append(ufloat(c, c_err))
        r2list.append(r2_value)

    latter_clist = clist[-int(len(clist)/2):]
    measured_gamma = sum(latter_clist)/len(latter_clist)
    # print("measured_gamma = ", measured_gamma)
    gamma_means = [clist[i].n for i in range(len(xvalues))]
    gamma_errors = [clist[i].s for i in range(len(xvalues))]
    ax1.errorbar(xvalues, gamma_means, yerr=gamma_errors, ls='none', capsize=3, color='k')

    if "nu_1_2" in file:
        ax1.axhline(-np.log(np.sqrt(2)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_1_3" in file or "nu_2_3" in file or "nu_1_3" in file:
        ax1.axhline(-np.log(np.sqrt(3)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_2_5" in file:
        ax1.axhline(-np.log(np.sqrt(5)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_4" in file or "nu_1_4" in file:
        ax1.axhline(-np.log(np.sqrt(4)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_7" in file:
        ax1.axhline(-np.log(np.sqrt(7)), color='r', linewidth=1, label="theory", zorder=2)

    ax1.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax1.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax1.set_xlabel("$\left(L_y/l_\mathrm{B}\\right)_\mathrm{min}$", fontsize=11)
    ax1.set_ylabel("$-\gamma_{\geq}$", fontsize=11)
    ax1.set_xlim([0, None])

    ####################################################################################################################

    ax2 = plt.subplot(upper_inner_grid[4], sharex=ax)

    ax2.scatter(xvalues, r2list, color='k', marker='x', s=10)
    ax2.axhline(0.99, color='r', linewidth=1, linestyle='dashed', label="threshold", zorder=2)

    ax2.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax2.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax2.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax2.set_ylabel("$R_{\geq}^2$", fontsize=11)
    ax2.set_xlim([0, None])
    ax2.set_ylim([None, 1])

    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_1_3"

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

    ax3 = plt.subplot(upper_inner_grid[1])
    ax3.tick_params('x', direction='in', bottom=True)

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

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if not error_bars:
                    ax3.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                            label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax3.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                            markersize=6)
                else:  # error_bars
                    ax3.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax3.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r',
                                marker=markers[flux_density_index], markersize=6)

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
            straight_line_of_best_fit(ax3, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    ax3.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

    ax3.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.53, 2), handletextpad=0, borderpad=0.4, framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax3.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax3.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax3.set_xlim([0, None])
    ax3.set_ylim(-1, 3)

    ####################################################################################################################

    ax4 = plt.subplot(upper_inner_grid[3], sharex=ax3)
    ax4.tick_params('x', direction='in', bottom=True)

    # plot the data by flux density
    LylB_total = []
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB_total.append(data_line[3])
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])

    clist = []
    r2list = []

    xvalues = sorted(LylB_total, key=float)[:len(LylB_total) - 2]

    # print(xvalues)

    for LylB_min in xvalues:
        # plot the line of best fit
        LylB = []
        SvN = []
        for i, val in enumerate(grouped_data):
            if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
                continue
            else:
                if val[3] >= LylB_min:
                    LylB.append(val[3])
                    SvN.append(val[4])
        m, m_err, c, c_err, r2_value = line_of_best_fit_values(LylB, SvN)
        clist.append(ufloat(c, c_err))
        r2list.append(r2_value)

    latter_clist = clist[-int(len(clist) / 2):]
    measured_gamma = sum(latter_clist) / len(latter_clist)
    # print("measured_gamma = ", measured_gamma)
    gamma_means = [clist[i].n for i in range(len(xvalues))]
    gamma_errors = [clist[i].s for i in range(len(xvalues))]
    ax4.errorbar(xvalues, gamma_means, yerr=gamma_errors, ls='none', capsize=3, color='k')

    if "nu_1_2" in file:
        ax4.axhline(-np.log(np.sqrt(2)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_1_3" in file or "nu_2_3" in file or "nu_1_3" in file:
        ax4.axhline(-np.log(np.sqrt(3)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_2_5" in file:
        ax4.axhline(-np.log(np.sqrt(5)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_4" in file or "nu_1_4" in file:
        ax4.axhline(-np.log(np.sqrt(4)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_7" in file:
        ax4.axhline(-np.log(np.sqrt(7)), color='r', linewidth=1, label="theory", zorder=2)

    ax4.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax4.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax4.set_xlabel("$\left(L_y/l_\mathrm{B}\\right)_\mathrm{min}$", fontsize=11)
    ax4.set_ylabel("$-\gamma_{\geq}$", fontsize=11)
    ax4.set_xlim([0, None])

    ####################################################################################################################

    ax5 = plt.subplot(upper_inner_grid[5], sharex=ax3)

    ax5.scatter(xvalues, r2list, color='k', marker='x', s=10)
    ax5.axhline(0.99, color='r', linewidth=1, linestyle='dashed', label="threshold", zorder=2)

    ax5.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax5.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax5.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax5.set_ylabel("$R_{\geq}^2$", fontsize=11)
    ax5.set_xlim([0, None])
    ax5.set_ylim([None, 1])

    ##################################################################################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_1_{filling}_total.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = True

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
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
        #                        LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
        LylB_outlier_values = [LylB_func(5 / 13, 7), LylB_func(2 / 9, 15), LylB_func(2 / 7, 15), LylB_func(5 / 14, 10),
                               LylB_func(1 / 4, 15), LylB_func(1 / 5, 15), LylB_func(3 / 8, 10), LylB_func(5 / 17, 8)]
    ####################################################################################################################

    ax6 = plt.subplot(middle_inner_grid[0])
    ax6.tick_params('x', direction='in', bottom=True)

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

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if not error_bars:
                    ax6.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                            label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax6.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                            markersize=6)
                else:  # error_bars
                    ax6.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax6.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k',
                                marker=markers[flux_density_index], markersize=6)

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

    threshold = 0.952

    for i in range(length):
        _, _, _, _, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > threshold:
            straight_line_of_best_fit(ax6, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), ypos=4.5)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    ax6.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

    ax6.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.53, 2.15), handletextpad=0, borderpad=0.4, framealpha=1,
              edgecolor='k', markerscale=1,
              fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax6.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax6.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax6.set_xlim([0, None])
    ax6.set_ylim([-1, None])

    ####################################################################################################################

    ax7 = plt.subplot(middle_inner_grid[2], sharex=ax6)
    ax7.tick_params('x', direction='in', bottom=True)

    # plot the data by flux density
    LylB_total = []
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB_total.append(data_line[3])
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])

    clist = []
    r2list = []

    xvalues = sorted(LylB_total, key=float)[:len(LylB_total) - 2]

    # print(xvalues)

    for LylB_min in xvalues:
        # plot the line of best fit
        LylB = []
        SvN = []
        for i, val in enumerate(grouped_data):
            if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
                continue
            else:
                if val[3] >= LylB_min:
                    LylB.append(val[3])
                    SvN.append(val[4])
        m, m_err, c, c_err, r2_value = line_of_best_fit_values(LylB, SvN)
        clist.append(ufloat(c, c_err))
        r2list.append(r2_value)

    latter_clist = clist[-int(len(clist) / 2):]
    measured_gamma = sum(latter_clist) / len(latter_clist)
    # print("measured_gamma = ", measured_gamma)
    gamma_means = [clist[i].n for i in range(len(xvalues))]
    gamma_errors = [clist[i].s for i in range(len(xvalues))]
    ax7.errorbar(xvalues, gamma_means, yerr=gamma_errors, ls='none', capsize=3, color='k')

    if "nu_1_2" in file:
        ax7.axhline(-np.log(np.sqrt(2)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_1_3" in file or "nu_2_3" in file or "nu_1_3" in file:
        ax7.axhline(-np.log(np.sqrt(3)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_2_5" in file:
        ax7.axhline(-np.log(np.sqrt(5)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_4" in file or "nu_1_4" in file:
        ax7.axhline(-np.log(np.sqrt(4)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_7" in file:
        ax7.axhline(-np.log(np.sqrt(7)), color='r', linewidth=1, label="theory", zorder=2)

    ax7.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax7.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax7.set_xlabel("$\left(L_y/l_\mathrm{B}\\right)_\mathrm{min}$", fontsize=11)
    ax7.set_ylabel("$-\gamma_{\geq}$", fontsize=11)
    ax7.set_xlim([0, None])

    ####################################################################################################################

    ax8 = plt.subplot(middle_inner_grid[4], sharex=ax6)

    ax8.scatter(xvalues, r2list, color='k', marker='x', s=10)
    ax8.axhline(threshold, color='r', linewidth=1, linestyle='dashed', label="threshold", zorder=2)

    ax8.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax8.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax8.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax8.set_ylabel("$R_{\geq}^2$", fontsize=11)
    ax8.set_xlim([0, None])
    ax8.set_ylim([None, 1])

    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"

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
        # LylB_outlier_values = [LylB_func(1 / 3, 4), LylB_func(1 / 4, 4), LylB_func(1 / 5, 4), LylB_func(1 / 6, 4),
        #                        LylB_func(2 / 7, 4), LylB_func(2 / 9, 4), LylB_func(2 / 13, 4), LylB_func(4 / 13, 4),
        #                        LylB_func(1 / 4, 6), LylB_func(1 / 5, 6), LylB_func(1 / 6, 6), LylB_func(2 / 7, 6),
        #                        LylB_func(2 / 9, 6), LylB_func(2 / 13, 6),
        #                        LylB_func(1 / 5, 8), LylB_func(1 / 6, 8), LylB_func(2 / 9, 8), LylB_func(2 / 13, 8)]
        # LylB_outlier_values = [LylB_func(1 / 7, 15), LylB_func(5 / 16, 9), LylB_func(2 / 7, 15), LylB_func(5 / 16, 8)]
        LylB_outlier_values = [LylB_func(5 / 13, 7), LylB_func(2 / 9, 15), LylB_func(2 / 7, 15), LylB_func(5 / 14, 10), LylB_func(1 / 4, 15), LylB_func(1 / 5, 15), LylB_func(3 / 8, 10), LylB_func(5 / 17, 8)]
        # LylB_outlier_values = [LylB_func(nphi, Ly) for nphi in [2/7] for Ly in range(20)] + [LylB_func(5 / 13, 7), LylB_func(2 / 9, 15), LylB_func(1 / 4, 15), LylB_func(5 / 17, 8)]
        print(LylB_outlier_values)

    ####################################################################################################################

    ax9 = plt.subplot(middle_inner_grid[1])
    ax9.tick_params('x', direction='in', bottom=True)

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

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if not error_bars:
                    ax9.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                             label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax9.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                             markersize=6)
                else:  # error_bars
                    ax9.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    # ax9.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='k',
                    #              marker=markers[flux_density_index], markersize=6)

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

    threshold = 0.99

    for i in range(length):
        _, _, _, _, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
        if r2value > threshold:
            straight_line_of_best_fit(ax9, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist(), ypos=0.75)
            critical_LylB = sorted_plot_data[0][0]
            break
        sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    ax9.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

    ax9.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.53, 2), handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax9.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax9.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax9.set_xlim([0, None])
    ax9.set_ylim([-1, None])

    ####################################################################################################################

    ax10 = plt.subplot(middle_inner_grid[3], sharex=ax9)
    ax10.tick_params('x', direction='in', bottom=True)

    # plot the data by flux density
    LylB_total = []
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB_total.append(data_line[3])
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])

    clist = []
    r2list = []

    xvalues = sorted(LylB_total, key=float)[:len(LylB_total) - 2]

    # print(xvalues)

    for LylB_min in xvalues:
        # plot the line of best fit
        LylB = []
        SvN = []
        for i, val in enumerate(grouped_data):
            if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
                continue
            else:
                if val[3] >= LylB_min:
                    LylB.append(val[3])
                    SvN.append(val[4])
        m, m_err, c, c_err, r2_value = line_of_best_fit_values(LylB, SvN)
        clist.append(ufloat(c, c_err))
        r2list.append(r2_value)

    latter_clist = clist[-int(len(clist) / 2):]
    measured_gamma = sum(latter_clist) / len(latter_clist)
    # print("measured_gamma = ", measured_gamma)
    gamma_means = [clist[i].n for i in range(len(xvalues))]
    gamma_errors = [clist[i].s for i in range(len(xvalues))]
    ax10.errorbar(xvalues, gamma_means, yerr=gamma_errors, ls='none', capsize=3, color='k')

    if "nu_1_2" in file:
        ax10.axhline(-np.log(np.sqrt(2)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_1_3" in file or "nu_2_3" in file or "nu_1_3" in file:
        ax10.axhline(-np.log(np.sqrt(3)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_2_5" in file:
        ax10.axhline(-np.log(np.sqrt(5)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_4" in file or "nu_1_4" in file:
        ax10.axhline(-np.log(np.sqrt(4)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_7" in file:
        ax10.axhline(-np.log(np.sqrt(7)), color='r', linewidth=1, label="theory", zorder=2)

    ax10.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax10.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax10.set_xlabel("$\left(L_y/l_\mathrm{B}\\right)_\mathrm{min}$", fontsize=11)
    ax10.set_ylabel("$-\gamma_{\geq}$", fontsize=11)
    ax10.set_xlim([0, None])

    ####################################################################################################################

    ax11 = plt.subplot(middle_inner_grid[5], sharex=ax9)

    ax11.scatter(xvalues, r2list, color='k', marker='x', s=10)
    ax11.axhline(threshold, color='r', linewidth=1, linestyle='dashed', label="threshold", zorder=2)

    ax11.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax11.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax11.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax11.set_ylabel("$R_{\geq}^2$", fontsize=11)
    ax11.set_xlim([0, None])
    ax11.set_ylim([None, 1])

    #######################################################################################################################################################################################

    model = "BosHofSqu1"
    filling = "nu_2_3"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_{filling}_total.out'

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

    ax12 = plt.subplot(lower_inner_grid[0])
    ax12.tick_params('x', direction='in', bottom=True)

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

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if not error_bars:
                    ax12.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                             label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax12.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                             markersize=6)
                else:  # error_bars
                    ax12.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax12.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r',
                                 marker=markers[flux_density_index], markersize=6)

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

    # for i in range(length):
    #     _, _, _, _, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
    #     if r2value > 0.99:
    #         straight_line_of_best_fit(ax12, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
    #         critical_LylB = sorted_plot_data[0][0]
    #         break
    #     sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    # ax12.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

    ax12.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.53, 2), handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax12.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax12.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax12.set_xlim([0, None])
    ax12.set_ylim([-1, None])

    ####################################################################################################################

    ax13 = plt.subplot(lower_inner_grid[2], sharex=ax12)
    ax13.tick_params('x', direction='in', bottom=True)

    # plot the data by flux density
    LylB_total = []
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB_total.append(data_line[3])
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])

    clist = []
    r2list = []

    xvalues = sorted(LylB_total, key=float)[:len(LylB_total) - 2]

    # print(xvalues)

    for LylB_min in xvalues:
        # plot the line of best fit
        LylB = []
        SvN = []
        for i, val in enumerate(grouped_data):
            if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
                continue
            else:
                if val[3] >= LylB_min:
                    LylB.append(val[3])
                    SvN.append(val[4])
        m, m_err, c, c_err, r2_value = line_of_best_fit_values(LylB, SvN)
        clist.append(ufloat(c, c_err))
        r2list.append(r2_value)

    latter_clist = clist[-int(len(clist) / 2):]
    measured_gamma = sum(latter_clist) / len(latter_clist)
    # print("measured_gamma = ", measured_gamma)
    gamma_means = [clist[i].n for i in range(len(xvalues))]
    gamma_errors = [clist[i].s for i in range(len(xvalues))]
    ax13.errorbar(xvalues, gamma_means, yerr=gamma_errors, ls='none', capsize=3, color='k')

    if "nu_1_2" in file:
        ax13.axhline(-np.log(np.sqrt(2)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_1_3" in file or "nu_2_3" in file or "nu_1_3" in file:
        ax13.axhline(-np.log(np.sqrt(3)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_2_5" in file:
        ax13.axhline(-np.log(np.sqrt(5)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_4" in file or "nu_1_4" in file:
        ax13.axhline(-np.log(np.sqrt(4)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_7" in file:
        ax13.axhline(-np.log(np.sqrt(7)), color='r', linewidth=1, label="theory", zorder=2)

    ax13.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    # ax13.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax13.set_xlabel("$\left(L_y/l_\mathrm{B}\\right)_\mathrm{min}$", fontsize=11)
    ax13.set_ylabel("$-\gamma_{\geq}$", fontsize=11)
    ax13.set_xlim([0, None])

    ####################################################################################################################

    ax14 = plt.subplot(lower_inner_grid[4], sharex=ax12)

    ax14.scatter(xvalues, r2list, color='k', marker='x', s=10)
    ax14.axhline(0.99, color='r', linewidth=1, linestyle='dashed', label="threshold", zorder=2)

    ax14.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    # ax14.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax14.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax14.set_ylabel("$R_{\geq}^2$", fontsize=11)
    ax14.set_xlim([0, None])
    ax14.set_ylim([None, 1])

    ####################################################################################################################

    model = "FerHofSqu1"
    filling = "nu_2_5"

    # specify the input file
    file = f'/home/bart/PycharmProjects/infinite_cylinder/logs/observables/{model}/{model}_Vrange_1_{filling}_total.out'

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

    ax15 = plt.subplot(lower_inner_grid[1])
    ax15.tick_params('x', direction='in', bottom=True)

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

    # plot the data by flux density
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                if not error_bars:
                    ax15.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                             label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax15.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                             markersize=6)
                else:  # error_bars
                    ax15.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                 label=f'${data_line[0]}/{data_line[1]}$', markersize=6)
                    ax15.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r',
                                 marker=markers[flux_density_index], markersize=6)

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

    # for i in range(length):
    #     _, _, _, _, r2value = line_of_best_fit(sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
    #     if r2value > 0.99:
    #         straight_line_of_best_fit(ax15, sorted_plot_data[:, 0].tolist(), sorted_plot_data[:, 1].tolist())
    #         critical_LylB = sorted_plot_data[0][0]
    #         break
    #     sorted_plot_data = np.delete(sorted_plot_data, 0, axis=0)

    ax15.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)

    ax15.legend(title='$n_\phi=p/q$', loc='upper center', bbox_to_anchor=(0.53, 2.15), handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)
    ax15.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax15.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax15.set_xlim([0, None])
    ax15.set_ylim([-1, None])

    ####################################################################################################################

    ax16 = plt.subplot(lower_inner_grid[3], sharex=ax15)
    ax16.tick_params('x', direction='in', bottom=True)

    # plot the data by flux density
    LylB_total = []
    for flux_density_index in range(len(flux_grouped_data)):
        LylB, SvN, SvN_error = [], [], []
        LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[
                2] < Ly_min:
                LylB_outliers.append(data_line[3])
                SvN_outliers.append(data_line[4])
                SvN_outliers_error.append(data_line[5])
            else:
                LylB_total.append(data_line[3])
                LylB.append(data_line[3])
                SvN.append(data_line[4])
                SvN_error.append(data_line[5])

    clist = []
    r2list = []

    xvalues = sorted(LylB_total, key=float)[:len(LylB_total) - 2]

    # print(xvalues)

    for LylB_min in xvalues:
        # plot the line of best fit
        LylB = []
        SvN = []
        for i, val in enumerate(grouped_data):
            if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
                continue
            else:
                if val[3] >= LylB_min:
                    LylB.append(val[3])
                    SvN.append(val[4])
        m, m_err, c, c_err, r2_value = line_of_best_fit_values(LylB, SvN)
        clist.append(ufloat(c, c_err))
        r2list.append(r2_value)

    latter_clist = clist[-int(len(clist) / 2):]
    measured_gamma = sum(latter_clist) / len(latter_clist)
    # print("measured_gamma = ", measured_gamma)
    gamma_means = [clist[i].n for i in range(len(xvalues))]
    gamma_errors = [clist[i].s for i in range(len(xvalues))]
    ax16.errorbar(xvalues, gamma_means, yerr=gamma_errors, ls='none', capsize=3, color='k')

    if "nu_1_2" in file:
        ax16.axhline(-np.log(np.sqrt(2)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_1_3" in file or "nu_2_3" in file or "nu_1_3" in file:
        ax16.axhline(-np.log(np.sqrt(3)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_2_5" in file:
        ax16.axhline(-np.log(np.sqrt(5)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_4" in file or "nu_1_4" in file:
        ax16.axhline(-np.log(np.sqrt(4)), color='r', linewidth=1, label="theory", zorder=2)
    elif "nu_3_7" in file:
        ax16.axhline(-np.log(np.sqrt(7)), color='r', linewidth=1, label="theory", zorder=2)

    ax16.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
                edgecolor='k', markerscale=1,
                fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax16.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax16.set_xlabel("$\left(L_y/l_\mathrm{B}\\right)_\mathrm{min}$", fontsize=11)
    ax16.set_ylabel("$-\gamma_{\geq}$", fontsize=11)
    ax16.set_xlim([0, None])

    ####################################################################################################################

    ax17 = plt.subplot(lower_inner_grid[5], sharex=ax15)

    ax17.scatter(xvalues, r2list, color='k', marker='x', s=10)
    ax17.axhline(0.99, color='r', linewidth=1, linestyle='dashed', label="threshold", zorder=2)

    ax17.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
                edgecolor='k', markerscale=1,
                fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    ax17.axvline(critical_LylB, color='g', ls='dashed', linewidth=1, zorder=2)
    ax17.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax17.set_ylabel("$R_{\geq}^2$", fontsize=11)
    ax17.set_xlim([0, None])
    ax17.set_ylim([None, 1])

    ##############################################################################################################################################################################

    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)

    plt.setp(ax6.get_xticklabels(), visible=False)
    plt.setp(ax7.get_xticklabels(), visible=False)
    plt.setp(ax9.get_xticklabels(), visible=False)
    plt.setp(ax10.get_xticklabels(), visible=False)

    plt.setp(ax12.get_xticklabels(), visible=False)
    plt.setp(ax13.get_xticklabels(), visible=False)
    plt.setp(ax15.get_xticklabels(), visible=False)
    plt.setp(ax16.get_xticklabels(), visible=False)

    ax.tick_params(axis='both', which='major', labelsize=10)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax3.tick_params(axis='both', which='major', labelsize=10)
    ax4.tick_params(axis='both', which='major', labelsize=10)
    ax5.tick_params(axis='both', which='major', labelsize=10)
    ax6.tick_params(axis='both', which='major', labelsize=10)
    ax7.tick_params(axis='both', which='major', labelsize=10)
    ax8.tick_params(axis='both', which='major', labelsize=10)
    ax9.tick_params(axis='both', which='major', labelsize=10)
    ax10.tick_params(axis='both', which='major', labelsize=10)
    ax11.tick_params(axis='both', which='major', labelsize=10)
    ax12.tick_params(axis='both', which='major', labelsize=10)
    ax13.tick_params(axis='both', which='major', labelsize=10)
    ax14.tick_params(axis='both', which='major', labelsize=10)
    ax15.tick_params(axis='both', which='major', labelsize=10)
    ax16.tick_params(axis='both', which='major', labelsize=10)
    ax17.tick_params(axis='both', which='major', labelsize=10)

    #fig.text(0.275, 0.96, "bosons", fontsize=16)
    #fig.text(0.7, 0.96, "fermions", fontsize=16)

    fig.text(0.075, 0.91, "(a) $\\nu=1/2$\n    error$\leq0.1\%$", fontsize=12)
    fig.text(0.5, 0.91, "(b) $\\nu=1/3$\n    error$\leq0.1\%$", fontsize=12)

    fig.text(0.075, 0.62, "(c) $\\nu=2/5$\n    no outliers", fontsize=12)
    fig.text(0.5, 0.62, "(d) $\\nu=2/5$\n    error$\leq0.1\%$", fontsize=12)

    fig.text(0.075, 0.33, "(e) $\\nu=2/3$\n    all data", fontsize=12)
    fig.text(0.5, 0.33, "(f) $\\nu=2/5$\n    all data", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TEE/figures/complete_data.png", bbox_inches='tight', dpi=300)
    plt.show()

