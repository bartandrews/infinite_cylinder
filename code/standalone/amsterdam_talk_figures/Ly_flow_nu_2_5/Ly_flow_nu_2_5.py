import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
from matplotlib.patches import ConnectionPatch
from numpy.polynomial.polynomial import polyfit
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker

if __name__ == '__main__':

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[0.75, 1])

    ####################################################################################################################

    ax1 = plt.subplot(gs[0, :])

    gs.update(hspace=0.3, wspace=0.28)

    Ly = []
    SvN = []
    Sinf = []

    with open('ent_scal.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            Ly.append(float(row[0]))
            SvN.append(float(row[1]))
            Sinf.append(float(row[2]))

    # ax1.plot(Ly[0], SvN[0], '.', marker='X', c='C8', label='$S_{\mathrm{vN}}\;(n_\phi=1/3)$', markersize=5, zorder=2, fillstyle='none')
    # ax1.plot(Ly[1], SvN[1], '.', marker='p', c='C8', label='$S_{\mathrm{vN}}\;(n_\phi=1/4)$', markersize=5, zorder=2, fillstyle='none')
    # ax1.plot(Ly[2], SvN[2], '.', marker='*', c='C8', label='$S_{\mathrm{vN}}\;(n_\phi=1/5)$', markersize=5, zorder=2, fillstyle='none')
    # ax1.plot(Ly[3], SvN[3], '.', marker='X', c='C8', markersize=5, zorder=2, fillstyle='none')
    # ax1.plot(Ly[4], SvN[4], '.', marker='p', c='C8', markersize=5, zorder=2, fillstyle='none')
    # ax1.plot(Ly[5], SvN[5], '.', marker='*', c='C8', markersize=5, zorder=2, fillstyle='none')
    # y = mx + c
    # parameters, V = np.polyfit(Ly, SvN, 1, cov=True)
    # m, c = parameters[0], parameters[1]
    # xvalues = np.arange(max(Ly) + 20)
    # ax1.plot(xvalues, m * xvalues + c, '-', c='C8', zorder=2)
    # print("SvN error in (m, c) = (",V[0][0], ",", V[1][1],")")
    # ax1.text(0.2, 2.2, "$S_\mathrm{{vN}}={gradient:.2f}L_y{intercept:.2f}$".format(gradient=m, intercept=c))

    ax1.plot(Ly[0], SvN[0], '.', marker='X', c='C8', label='$n_\phi=1/3$', markersize=5, zorder=2, fillstyle='none')
    ax1.plot(Ly[1], SvN[1], '.', marker='p', c='C8', label='$n_\phi=1/4$', markersize=5, zorder=2, fillstyle='none')
    ax1.plot(Ly[2], SvN[2], '.', marker='*', c='C8', label='$n_\phi=1/5$', markersize=5, zorder=2, fillstyle='none')

    # y = mx + c
    # parameters, V = np.polyfit(Ly, Sinf, 1, cov=True)
    # m, c = parameters[0], parameters[1]
    m = 0.204987467
    c = -0.877999554
    xvalues = np.arange(max(Ly) + 20)
    ax1.plot(xvalues, m * xvalues + c, '-', c='C8', zorder=2)
    # print("Sinf error in (m, c) = (",V[0][0], ",", V[1][1],")")
    ax1.text(0.2, 0.5, "$S_{{\\mathrm{{vN}}}}={gradient:.2f}(L_y/l_\mathrm{{B}})-({intercept:.2f}\pm 0.08)$".format(gradient=m, intercept=abs(c)))


    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.6), handletextpad=0, borderpad=0.4, framealpha=1, edgecolor='k', markerscale=1,
               fontsize=10, ncol=3, labelspacing=0, columnspacing=0)

    ax1.set_xticks(np.arange(0, 16.1, 4))
    ax1.set_xlim([0, 12])
    ax1.set_xlabel("$L_y/l_\mathrm{B}$", fontsize=11)
    ax1.set_yticks(np.arange(-1, 1.1, 1))
    ax1.set_ylim([-1, 1])
    ax1.set_ylabel("$S_{\mathrm{vN}}$", fontsize=11)

    ax1.axhline(-0.805, color='k', linewidth=0.5, ls='--')
    #ax1.axvline(9.33072, color='k', linewidth=0.5, ls='-', zorder=1)
    ax1.axvline(7.77560150778107, color='k', linewidth=0.5, ls='-', zorder=1)

    ax1.tick_params(axis="x", labelsize=10)
    ax1.tick_params(axis="y", labelsize=10)

    fig.subplots_adjust(top=0.8)

    # ax1.xaxis.labelpad = -3
    ax1.xaxis.set_label_coords(0.5, -0.2)

    ax1.set_aspect(1)

    ####################################################################################################################

    ax2 = plt.subplot(gs[1, 0])

    x = []
    y = []
    z = []

    with open('ent_spec_mom_Ly_5.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(3, -3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax2.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='{}'.format(value))
        else:
            ax2.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='{}'.format(value))

    ax2.set_yticks(np.arange(0, 15.1, 5))
    ax2.set_ylim([0, 15])
    ax2.legend(loc='center left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2,
                    fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))
    ax2.set_xlim([-np.pi/3, np.pi/3])
    ax2.set_xlabel("$k_\\alpha / \pi$", fontsize=11)
    ax2.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ax2.tick_params(axis="x", labelsize=10)
    ax2.tick_params(axis="y", labelsize=10)

    ax2.plot([-0.1, 1*np.pi/3], [0, 10], color='k', linestyle=':', linewidth=1)
    ax2.plot([0.65, 1*np.pi/3], [0, 3.5], color='k', linestyle=':', linewidth=1)

    # ax2.axhline(1, color='k', linewidth=0.5, ls='--')
    # ax2.axvline(2, color='k', linewidth=0.5, ls='--')
    # ax2.axhline(0.153205653777937, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(2.073942791479203, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(4.410048930142911, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(5.700361870617707, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(7.870452206203998, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(8.400756721645187, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(8.802885632428440, color='k', linewidth=0.5, ls=':')

    fig.text(0.14, 0.4625, '$L_y=5$ sites', fontsize=11, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    # gs.update(wspace=0.25, hspace=0.2)

    ####################################################################################################################

    # ax3 = plt.subplot(gs[1, 1])
    #
    # x = []
    # y = []
    # z = []
    #
    # with open('ent_spec_mom_Ly_10.dat', 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     for row in plots:
    #         x.append(float(row[1]))
    #         y.append(float(row[2]))
    #         z.append(int(row[0]))
    #
    # for value in np.linspace(3, -3, 7, dtype=int):
    #     xvalue = []
    #     yvalue = []
    #     for i in range(len(x)):
    #         if z[i] == value:
    #             xvalue.append(x[i])
    #             yvalue.append(y[i])
    #     if value!=0:
    #         ax3.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4)%10), label='{}'.format(value))
    #     else:
    #         ax3.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='{}'.format(value))
    #
    # ax3.set_yticks(np.arange(0, 15.1, 5))
    # ax3.set_ylim([0, 15])
    # # ax3.legend(loc='center left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2,
    # #             fontsize=12, ncol=1, bbox_to_anchor=(1, 1))
    # ax3.set_xlim([-np.pi/3, np.pi/3])
    # ax3.set_xlabel("$k_\\alpha / \pi$", fontsize=11)
    # # ax3.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)
    #
    # ax3.tick_params(axis="x", labelsize=10)
    # ax3.tick_params(axis="y", labelsize=10)
    #
    # ax3.plot([-0.25, np.pi/3], [0, 13], color='k', linestyle=':', linewidth=1)
    # ax3.plot([0.65, np.pi / 3], [0, 4], color='k', linestyle=':', linewidth=1)
    #
    # # ax3.axvline(1, color='k', linewidth=0.5, ls='--')
    # # ax3.axvline(2, color='k', linewidth=0.5, ls='--')
    #
    # fig.text(0.574, 0.4625, '$L_y=10$ sites', fontsize=11, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=1))
    #
    # # gs.update(wspace=0.25, hspace=0)
    # plt.setp(ax3.get_yticklabels(), visible=False)
    # ax3.tick_params(axis='y', which='both', length=0)

    ####################################################################################################################

    Ly6_con = ConnectionPatch(xyA=(7.77560150778107, -1), xyB=(1, 15), coordsA="data", coordsB="data",
                              axesA=ax1, axesB=ax2, connectionstyle="angle3,angleA=-90,angleB=10", arrowstyle='->',
                              facecolor='k', edgecolor='k')
    # Ly9_con = ConnectionPatch(xyA=(13.99608, -1), xyB=(0.45, 15), coordsA="data", coordsB="data",
    #                           axesA=ax1, axesB=ax3, connectionstyle="angle3,angleA=30,angleB=90", arrowstyle='->', facecolor='k', edgecolor='k')

    ax1.add_artist(Ly6_con)
    # ax1.add_artist(Ly9_con)

    # fig.text(0.035, 0.8, "(a)", color="k", fontsize=12)
    # fig.text(0.035, 0.47, "(b)", color="k", fontsize=12)
    # fig.text(0.48, 0.49, "(c)", color="k", fontsize=12)

    fig.text(0.035, 0.82, "$l_\mathrm{B}\sim n_\phi^{-1/2}$", color="k", fontsize=14)
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    fig.text(0.11, 0.58, "o", color="k", fontsize=25)

    plt.savefig("/home/bart/Documents/presentations/2021_03_18/figures/Ly_flow_nu_2_5_new.png", bbox_inches='tight', dpi=300)
    plt.show()