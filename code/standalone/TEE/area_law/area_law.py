import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import groupby
# import sys
import math
from scipy import stats
from uncertainties import unumpy, ufloat


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


def line_of_best_fit_values(LylB_list, SvN_list):

    parameters, cov = np.polyfit(LylB_list, SvN_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(LylB_list, SvN_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    return m, m_err, c, c_err, r2_value


if __name__ == '__main__':

    # specify the input file
    file = '/home/bart/PycharmProjects/infinite_cylinder/logs/observables/BosHofSqu1/BosHofSqu1_nu_2_3_total.out'

    # plot with error bars?
    error_bars = True

    # identify the outliers?
    identify_outliers = False

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

    fig = plt.figure()
    ax = plt.subplot(111)

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
                    ax.plot(LylB, SvN, '.', marker=markers[flux_density_index], label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=6)
                    ax.plot(LylB_outliers, SvN_outliers, '.', color='k', marker=markers[flux_density_index], markersize=6)
                else:  # error_bars
                    ax.errorbar(LylB, SvN, yerr=SvN_error, ls='none', capsize=3, marker=markers[flux_density_index], label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=6)
                    ax.errorbar(LylB_outliers, SvN_outliers, yerr=SvN_outliers_error, ls='none', capsize=3, color='r', marker=markers[flux_density_index], markersize=6)

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
            continue
        else:
            if val[3] > 10:
                LylB.append(val[3])
                SvN.append(val[4])

    line_of_best_fit(LylB, SvN)

    # ax.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
    #           edgecolor='k', markerscale=1,
    #           fontsize=10, ncol=3, labelspacing=0, columnspacing=0)
    ax.set_title(file)
    ax.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)
    ax.set_xlim(0)
    ax.set_ylim(-1)

    ####################################################################################################################
    #
    # #fig = plt.figure()
    # #ax2 = plt.subplot(111)
    #
    # # plot the data by flux density
    # LylB_total = []
    # for flux_density_index in range(len(flux_grouped_data)):
    #     LylB, SvN, SvN_error = [], [], []
    #     LylB_outliers, SvN_outliers, SvN_outliers_error = [], [], []
    #     for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
    #         if any(math.isclose(j, data_line[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or data_line[2] < Ly_min:
    #             LylB_outliers.append(data_line[3])
    #             SvN_outliers.append(data_line[4])
    #             SvN_outliers_error.append(data_line[5])
    #         else:
    #             LylB_total.append(data_line[3])
    #             LylB.append(data_line[3])
    #             SvN.append(data_line[4])
    #             SvN_error.append(data_line[5])
    #
    # clist = []
    #
    # xvalues = sorted(LylB_total, key=float)[:len(LylB_total)-2]
    #
    # print(xvalues)
    #
    # for LylB_min in xvalues:
    #     # plot the line of best fit
    #     LylB = []
    #     SvN = []
    #     for i, val in enumerate(grouped_data):
    #         if any(math.isclose(j, val[3], rel_tol=1e-5) is True for j in LylB_outlier_values) or val[2] < Ly_min:
    #             continue
    #         else:
    #             if val[3] >= LylB_min:
    #                 LylB.append(val[3])
    #                 SvN.append(val[4])
    #     m, m_err, c, c_err, r2_value = line_of_best_fit_values(LylB, SvN)
    #     clist.append(ufloat(c, c_err))
    #
    # latter_clist = clist[-int(len(clist)/2):]
    # measured_gamma = sum(latter_clist)/len(latter_clist)
    # print("measured_gamma = ", measured_gamma)
    # gamma_means = [clist[i].n for i in range(len(xvalues))]
    # gamma_errors = [clist[i].s for i in range(len(xvalues))]
    # ax.errorbar(xvalues, gamma_means, yerr=gamma_errors, ls='none', capsize=3)
    #
    # ax.axvline(xvalues[int(len(clist) / 2)], color='k', linewidth=0.5, ls='--')
    # # ax.axhline(-np.log(np.sqrt(2)), color='r', linewidth=2, ls=(0, (5, 10)), label="theory", zorder=2)
    # ax.axhline(measured_gamma.n, color='k', linewidth=1, ls='-', label="measured", zorder=1)
    #
    # ax.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
    #           edgecolor='k', markerscale=1,
    #           fontsize=10, ncol=4, labelspacing=0, columnspacing=0)
    #
    # fig.text(0.6, 0.4, "$\\bar{{\gamma}}_{{>}}={g:.3f}\pm{g_err:.3f}$".format(g=abs(measured_gamma.n), g_err=measured_gamma.s))

    plt.show()

