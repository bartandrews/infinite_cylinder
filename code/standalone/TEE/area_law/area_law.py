import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
from matplotlib.patches import ConnectionPatch
from numpy.polynomial.polynomial import polyfit
from mpl_toolkits.mplot3d import Axes3D
from itertools import groupby
import sys
from fractions import Fraction as Frac
import itertools

if __name__ == '__main__':

    file = '/home/bart/PycharmProjects/infinite_cylinder/logs/observables/BosHofSqu1/log_analysis.out'

    fig = plt.figure()
    ax = plt.subplot(111)

    ####################################################################################################################

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
            line_count += 1

    # group the list into lines of the file
    i = 0
    grouped_data = []
    while i < len(data):
        grouped_data.append(data[i:i + 5])
        i += 5

    # group data by flux density
    flux_grouped_data = [list(i) for j, i in groupby(grouped_data, key=lambda a: a[0]/a[1])]

    # plot the data by flux density
    marker = itertools.cycle(((3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0),
                              (4, 1, 0), (5, 1, 0), (6, 1, 0),
                              (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P'))
    for flux_density_index in range(len(flux_grouped_data)):
        LylB = []
        SvN = []
        for i, data_line in enumerate(flux_grouped_data[flux_density_index]):
            LylB.append(data_line[2])
            SvN.append(data_line[3])
            if i == len(flux_grouped_data[flux_density_index]) - 1:
                ax.plot(LylB, SvN, '.', marker=next(marker), label=f'$n_\phi={data_line[0]}/{data_line[1]}$', markersize=6)

    # plot the line of best fit
    LylB = []
    SvN = []
    for i, val in enumerate(grouped_data):
        LylB.append(grouped_data[i][2])
        SvN.append(grouped_data[i][3])

    parameters, cov = np.polyfit(LylB, SvN, 1, cov=True)
    m, c = parameters[0], parameters[1]
    print("SvN error in (m, c) = (", np.sqrt(cov[0][0]), ",", np.sqrt(cov[1][1]), ")")
    xvalues = np.arange(max(LylB))
    ax.plot(xvalues, m * xvalues + c, '-', c='k', zorder=0)
    ax.text(10, 0, "$S_\mathrm{{vN}}={gradient:.2f}(L_y/l_\mathrm{{B}})-({intercept:.2f}\pm {cerror:.2f})$".format(
                 gradient=m, intercept=abs(c), cerror=np.sqrt(cov[1][1])))

    ax.legend(loc='upper left', handletextpad=0, borderpad=0.4, framealpha=1,
               edgecolor='k', markerscale=1,
               fontsize=10, ncol=3, labelspacing=0, columnspacing=0)

    ax.set_title(file)
    ax.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax.set_ylabel("$S_\mathrm{{vN}}$", fontsize=11)

    ax.set_xlim(0)

    ####################################################################################################################

    plt.show()
