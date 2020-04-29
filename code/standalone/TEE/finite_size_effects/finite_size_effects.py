import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import groupby
# import sys
import math
from scipy import stats
import matplotlib.gridspec as gridspec


##################
# Ly/lB function #
##################


def LylB_func(nphi, Ly):
    return np.sqrt(2*np.pi*nphi)*Ly


#############################
# line of best fit function #
#############################


def line_of_best_fit(axis, LylB_list, SvN_list):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, c = parameters[0], parameters[1]
    m_err = np.sqrt(cov[0][0])
    c_err = np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    print("SvN = m*(Ly/lB) + c")
    print(f"(m, m_err, c, c_err) = ({m:.3f}, {m_err:.3f}, {c:.3f}, {c_err:.3f})")
    xvalues = np.arange(max(LylB_list) + 1)
    axis.plot(xvalues, m * xvalues + c, '-', c='k', zorder=0)
    # fig.text(0.4, 0.2, "$S_\mathrm{{vN}}={gradient:.3f}(L_y/l_\mathrm{{B}})+({intercept:.3f}\pm {cerror:.3f})$\n$R^2={rsquared:.10f}$".format(
    #     gradient=m, intercept=c, cerror=np.sqrt(cov[1][1]), rsquared=r2_value))

    return m, m_err, c, c_err, r2_value


if __name__ == '__main__':

    marker_size = 6

    ####################################################################################################################

    fig = plt.figure(figsize=(13.75, 2.5))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1], wspace=0.5)

    ####################################################################################################################

    ax1 = plt.subplot(gs[0])

    # specify the input file
    file = '/home/bart/PycharmProjects/infinite_cylinder/logs/observables/BosHofSqu1/BosHofSqu1_nu_1_2_accepted.out'

    # plot with error bars?
    error_bars = False

    # identify the outliers?
    identify_outliers = False

    # set Ly_min
    Ly_min = 4

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

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(float(row[2]))
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
                    ax1.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                             label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=marker_size)
                    ax1.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                             markersize=marker_size)
                else:  # error_bars
                    ax1.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                 label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=marker_size)
                    ax1.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r',
                                 marker=markers[flux_density_index], markersize=marker_size)

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    m, m_err, c, c_err, r2 = line_of_best_fit(ax1, LylB, SvN)

    ax1.axhline(c, color='k', linewidth=0.5, ls='--')
    ax1.text(2.5, c + 0.1, f"$\gamma={abs(c):.3f}\pm{c_err:.3f}$")
    ax1.text(1, 1.9, f"$R^2={r2:.3f}$")

    ax1.legend(loc='center', bbox_to_anchor=(2.75, -0.45), handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=8, labelspacing=0, columnspacing=0)
    ax1.set_title("$L_y\geq 4$")
    ax1.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax1.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax1.set_xlim(0)

    ####################################################################################################################

    ax2 = plt.subplot(gs[1])

    # specify the input file
    file = '/home/bart/PycharmProjects/infinite_cylinder/logs/observables/BosHofSqu1/BosHofSqu1_nu_1_2_accepted.out'

    # plot with error bars?
    error_bars = False

    # identify the outliers?
    identify_outliers = False

    # set Ly_min
    Ly_min = 6

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

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(float(row[2]))
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
                    ax2.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                             label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=marker_size)
                    ax2.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                             markersize=marker_size)
                else:  # error_bars
                    ax2.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                 label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=marker_size)
                    ax2.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r',
                                 marker=markers[flux_density_index], markersize=marker_size)

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    m, m_err, c, c_err, r2 = line_of_best_fit(ax2, LylB, SvN)

    ax2.axhline(c, color='k', linewidth=0.5, ls='--')
    ax2.text(2.5, c + 0.1, f"$\gamma={abs(c):.3f}\pm{c_err:.3f}$")
    ax2.text(1, 1.9, f"$R^2={r2:.3f}$")

    # ax2.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
    #           edgecolor='k', markerscale=1,
    #           fontsize=10, ncol=3, labelspacing=0, columnspacing=0)
    ax2.set_title("$L_y\geq 6$")
    ax2.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax2.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax2.set_xlim(0)

    ####################################################################################################################

    ax3 = plt.subplot(gs[2])

    # specify the input file
    file = '/home/bart/PycharmProjects/infinite_cylinder/logs/observables/BosHofSqu1/BosHofSqu1_nu_1_2_accepted.out'

    # plot with error bars?
    error_bars = False

    # identify the outliers?
    identify_outliers = False

    # set Ly_min
    Ly_min = 8

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

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(float(row[2]))
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
                             label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=marker_size)
                    ax3.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                             markersize=marker_size)
                else:  # error_bars
                    ax3.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                 label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=marker_size)
                    ax3.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r',
                                 marker=markers[flux_density_index], markersize=marker_size)

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    m, m_err, c, c_err, r2 = line_of_best_fit(ax3, LylB, SvN)

    ax3.axhline(c, color='k', linewidth=0.5, ls='--')
    ax3.text(2.5, c + 0.1, f"$\gamma={abs(c):.3f}\pm{c_err:.3f}$")
    ax3.text(1, 1.9, f"$R^2={r2:.3f}$")

    # ax3.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
    #           edgecolor='k', markerscale=1,
    #           fontsize=10, ncol=3, labelspacing=0, columnspacing=0)
    ax3.set_title("$L_y\geq 8$")
    ax3.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax3.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax3.set_xlim(0)

    ####################################################################################################################

    ax4 = plt.subplot(gs[3])

    # specify the input file
    file = '/home/bart/PycharmProjects/infinite_cylinder/logs/observables/BosHofSqu1/BosHofSqu1_nu_1_2_accepted.out'

    # plot with error bars?
    error_bars = False

    # identify the outliers?
    identify_outliers = False

    # set Ly_min
    Ly_min = 10

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

    # append data from file to a list
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        line_count = 0
        for row in plots:
            data.append(int(row[0]))
            data.append(int(row[1]))
            data.append(float(row[2]))
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
                    ax4.plot(LylB, SvN, '.', marker=markers[flux_density_index],
                            label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=marker_size)
                    ax4.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index],
                            markersize=marker_size)
                else:  # error_bars
                    ax4.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index],
                                label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=marker_size)
                    ax4.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r',
                                marker=markers[flux_density_index], markersize=marker_size)

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            LylB.append(val[3])
            SvN.append(val[4])

    m, m_err, c, c_err, r2 = line_of_best_fit(ax4, LylB, SvN)

    ax4.axhline(c, color='k', linewidth=0.5, ls='--')
    ax4.text(2.5, c+0.1, f"$\gamma={abs(c):.3f}\pm{c_err:.3f}$")
    ax4.text(1, 1.9, f"$R^2={r2:.3f}$")

    # ax4.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
    #           edgecolor='k', markerscale=1,
    #           fontsize=10, ncol=3, labelspacing=0, columnspacing=0)
    ax4.set_title("$L_y\geq 10$")
    ax4.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax4.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax4.set_xlim(0)
    
    ####################################################################################################################

    fig.text(0.08, 0.9, "(a)", fontsize=12)
    fig.text(0.29, 0.9, "(b)", fontsize=12)
    fig.text(0.5, 0.9, "(c)", fontsize=12)
    fig.text(0.71, 0.9, "(d)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TEE/figures/finite_size_effects.png", bbox_inches='tight', dpi=300)
    plt.show()
