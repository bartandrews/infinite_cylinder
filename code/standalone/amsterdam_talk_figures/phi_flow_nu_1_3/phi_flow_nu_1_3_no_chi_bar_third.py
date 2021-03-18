import numpy as np
import matplotlib
from matplotlib.colorbar import Colorbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
from matplotlib.patches import ConnectionPatch
from numpy.polynomial.polynomial import polyfit
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

    ####################################################################################################################

    ax1 = plt.subplot(gs[:, 0], anchor=(0, 0.85))

    gs.update(hspace=0)

    phi_IQH = []
    IQH_charge = []
    phi_FQH = []

    # an array of parameters, each of our curves depend on a specific
    # value of parameters
    chi_range = np.linspace(50, 150, 3)

    # norm is a class which, when called, can normalize data into the
    # [0.0, 1.0] interval.
    norm = matplotlib.colors.Normalize(
        vmin=np.min(chi_range),
        vmax=np.max(chi_range))

    # choose a colormap
    # c_m = matplotlib.cm.cool
    # easily visible colors = ["red", "blue"]
    # presentation colors = ["white", "C9"]
    c_m = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red", "blue"])

    # create a ScalarMappable and initialize a data structure
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    kappa_range = np.linspace(0, 2, 3)
    FQH_charge = {}
    for kappa in kappa_range:
        FQH_charge[kappa] = []

    FQH_charge_kappa_limit = {}
    for chi in chi_range:
        FQH_charge_kappa_limit[chi] = []

    with open('charge_pump_IQH.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            phi_IQH.append(float(row[0]))
            IQH_charge.append(float(row[1]))

    with open('charge_pump_FQH_kappa.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            phi_FQH.append(float(row[0]))
            for i, kappa in enumerate(kappa_range):
                FQH_charge[kappa].append(float(row[i+1]))

    with open('charge_pump_FQH_kappa_0.3.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            # phi_FQH.append(float(row[0]))
            for i, chi in enumerate(chi_range):
                FQH_charge_kappa_limit[chi].append(float(row[i+1]))

    FQH_charge_kappa_04 = []
    with open('charge_pump_FQH_kappa_0.4.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            FQH_charge_kappa_04.append(float(row[1]))

    ax1.plot(phi_IQH, IQH_charge, marker='o', linestyle='-', c='C8', label='$\\mathbf{IQH}$', markersize=5, mec='k',
             zorder=1)

    ax1.plot(phi_FQH, FQH_charge[0], marker='d', linestyle='-', color='C5', label='FQH ($\kappa=0)$', markersize=5,
             zorder=1)
    ax1.plot(phi_FQH, FQH_charge[1], marker='>', linestyle='-', color='C1', label='FQH ($\kappa=0.1)$', markersize=5,
             zorder=2)
    ax1.plot(phi_FQH, FQH_charge[2], marker='<', linestyle='-', color='C2', label='FQH ($\kappa=0.2)$', markersize=5,
             zorder=3)

    for i, chi in enumerate(chi_range):
        if i == len(chi_range)-1:
            ax1.plot(phi_FQH, FQH_charge_kappa_limit[chi], marker='s', linestyle='-', color="blue", label='$\\mathbf{FQH}$ ($\kappa=0.3$)', markersize=5, zorder=5)

    ax1.plot(phi_FQH, FQH_charge_kappa_04, marker='^', linestyle='-', color='C6', label='FQH ($\kappa=0.4)$', markersize=5,
             zorder=3)

    # cbax = plt.subplot()  # Place it where it should be.
    # # --------------------------------------------------------
    # cb = Colorbar(mappable=s_m, ax=cbax, orientation='horizontal', ticklocation='top')
    # cb.set_label(r'Colorbar !', labelpad=10)

    # cbax = fig.add_axes([0.125, 0.75+0.12, 0.343, 0.025])
    # cb = plt.colorbar(s_m, orientation='horizontal', cax=cbax, ticklocation='top', ticks=[50, 100, 150])
    # cb.set_label("$\chi$ for FQH ($\kappa=0.3$)", fontsize=11, labelpad=7)

    # ax1.legend(bbox_to_anchor = [0.62, 0.53], handletextpad=0.2, borderpad=0.4, framealpha=1, edgecolor='k', markerscale=1, fontsize=10)
    ax1.legend(bbox_to_anchor=(0.45, -0.61), loc="lower center", handletextpad=0.2, columnspacing=0.25, borderpad=0.4, framealpha=1, edgecolor='k',
               markerscale=1, fontsize=10, ncol=2)

    ax1.set_xticks(np.arange(0, 3.1, 1))
    ax1.set_xlim([0, 3])
    ax1.set_xlabel("$\Phi_x / 2\pi$", fontsize=11)
    ax1.set_yticks(np.arange(-3, 0.1, 1))
    ax1.set_ylim([-3, 0.1])
    ax1.set_ylabel("$\langle Q_{\mathrm{L}}(\Phi_x) \\rangle$", fontsize=11)

    ax1.axhline(-1, color='k', linewidth=0.5, ls='--')
    ax1.axhline(-2, color='k', linewidth=0.5, ls='--')
    ax1.axvline(1, color='k', linewidth=0.5, ls='--')
    ax1.axvline(2, color='k', linewidth=0.5, ls='--')

    ax1.tick_params(axis="x", labelsize=10)
    ax1.tick_params(axis="y", labelsize=10)

    ax1.set_aspect(1)

    ############################################

    left, bottom, width, height = [0.11, 0.215+0.12, 0.25, 0.25]  # left = 0.12
    ax5 = fig.add_axes([left, bottom, width, height], projection='3d')

    r = 5

    end = 6
    point_size = 2

    for y in range(0, end, 2):
        theta = np.linspace(0, 2 * np.pi, 11)
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='b', s=point_size)

        theta = np.linspace((1 / 20) * 2 * np.pi, (1 + 1 / 20) * 2 * np.pi, 11)
        y = y + 1 / 3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='r', s=point_size)

        theta = np.linspace((1 / 20) * 2 * np.pi, (1 + 1 / 20) * 2 * np.pi, 11)
        y = y + 2 / 3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='b', s=point_size)

        theta = np.linspace(0, 2 * np.pi, 11)
        y = y + 1 / 3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='r', s=point_size)

    # cut
    xx, zz = np.meshgrid(np.linspace(-7, 7, 10), np.linspace(-7, 7, 10))
    yy = end/2
    ax5.plot_surface(xx, yy, zz, alpha=0.25)

    # arrow
    ax5.quiver(0, end, 0, 0, -(end + end/4), 0, arrow_length_ratio=0.06, color='k', lw=0.5)
    ax5.text2D(-0.078, 0.006, '$\Phi_x$', fontsize=11)

    # cylinder
    angle = np.linspace(0, 2 * np.pi, 21)
    length = np.linspace(0, end-2/3, 13)
    theta, y = np.meshgrid(angle, length)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax5.plot_surface(x, y, z, alpha=0.25)

    # ring
    theta = np.linspace(0, 2 * np.pi, 21)
    y = np.repeat(end, 21)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax5.plot(x, y, z, label='parametric curve', c='k', lw=0.5)

    theta = np.linspace(0, 2 * np.pi, 1)
    y = np.repeat(end, 1)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax5.scatter(x, y, z, label='parametric curve', c='k', marker='v', s=5)
    ax5.text2D(0.05, 0.03, '$L_y$', fontsize=11)

    ax5.view_init(elev=10., azim=-15)
    ax5.axis('off')
    ax5.patch.set_alpha(0)

    ####################################################################################################################

    ax2 = plt.subplot(gs[0, 1])

    x = []
    y = []
    z = []

    with open('ent_spec_flow_IQH.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(4, -4, 9, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        ax2.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value+4)%10), label='{}'.format(value))

    ax2.set_yticks(np.arange(0, 8.1, 2))
    ax2.set_ylim([0, 8])
    ax2.legend(loc='center left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2,
                    fontsize=12, ncol=1, bbox_to_anchor=(1, 0))
    ax2.set_xlim([0, 3])
    # ax2.set_xlabel("$\Phi_y / 2\pi$", fontsize=11)
    ax2.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ax2.tick_params(axis="x", labelsize=10)
    ax2.tick_params(axis="y", labelsize=10)

    ax2.axvline(1, color='k', linewidth=0.5, ls='--')
    ax2.axvline(2, color='k', linewidth=0.5, ls='--')
    # ax2.axhline(0.153205653777937, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(2.073942791479203, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(4.410048930142911, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(5.700361870617707, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(7.870452206203998, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(8.400756721645187, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(8.802885632428440, color='k', linewidth=0.5, ls=':')

    ax2.arrow(0, 0, 0, 8, width=0.1, head_width=0, fill=False, zorder=3)
    ax2.arrow(1, 0, 0, 8, width=0.1, head_width=0, fill=False, zorder=3)
    ax2.arrow(2, 0, 0, 8, width=0.1, head_width=0, fill=False, zorder=3)
    ax2.arrow(3, 0, 0, 8, width=0.1, head_width=0, fill=False, zorder=3)

    ax2.tick_params('y', which="major", direction='in', right=True)

    fig.text(0.57, 0.86, 'IQH', fontsize=11, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    gs.update(wspace=0.25, hspace=0.2)

    ####################################################################################################################

    ax3 = plt.subplot(gs[1, 1])

    x = []
    y = []
    z = []

    with open('ent_spec_flow_FQH.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(4, -4, 9, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        ax3.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4)%10), label='{}'.format(value))

    ax3.set_yticks(np.arange(0, 8, 2))
    ax3.set_ylim([0, 8])
    # ax3.legend(loc='center left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2,
    #             fontsize=12, ncol=1, bbox_to_anchor=(1, 1))
    ax3.set_xlim([0, 3])
    ax3.set_xlabel("$\Phi_x / 2\pi$", fontsize=11)
    ax3.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ax3.tick_params(axis="x", labelsize=10)
    ax3.tick_params(axis="y", labelsize=10)

    ax3.axvline(1, color='k', linewidth=0.5, ls='--')
    ax3.axvline(2, color='k', linewidth=0.5, ls='--')

    ax3.arrow(0, 0, 0, 8, width=0.1, head_width=0, fill=False, zorder=3)
    ax3.arrow(3, 0, 0, 8, width=0.1, head_width=0, fill=False, zorder=3)

    ax3.tick_params('y', which="major", direction='in', right=True)

    fig.text(0.57, 0.475, 'FQH at $\kappa=0.3$', fontsize=11, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    gs.update(wspace=0.25, hspace=0)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ####################################################################################################################

    # IQH_con = ConnectionPatch(xyA=(2, -2.1), xyB=(0, 1), coordsA="data", coordsB="data",
    #                       axesA=ax1, axesB=ax3, connectionstyle="angle3,angleA=-90,angleB=180", arrowstyle='fancy', facecolor='C8', edgecolor='C8')
    # FQH_con = ConnectionPatch(xyA=(2, -0.4), xyB=(0, 5), coordsA="data", coordsB="data",
    #                          axesA=ax1, axesB=ax2, connectionstyle="angle3,angleA=90,angleB=180", arrowstyle='fancy',
    #                          facecolor='C9', edgecolor='C9')
    # ax1.add_artist(IQH_con)
    # ax1.add_artist(FQH_con)

    # fig.text(0.035, 0.83, "(a)", color="k", fontsize=12)
    # fig.text(0.49, 0.875, "(b)", color="k", fontsize=12)
    # fig.text(0.49, 0.49, "(c)", color="k", fontsize=12)

    fig.text(0.115, 0.87, "Hall response: $\sigma_\mathrm{H}=\\frac{e^2}{h}C\\nu$", color="k", fontsize=14)

    fig.text(0.34, 0.74, "$C\\nu=1/3$", color="k", fontsize=12)
    fig.text(0.145, 0.63, "$C\\nu=1$", color="k", fontsize=12)

    fig.text(0.45, 0.65, "o", color="k", fontsize=25)
    fig.text(0.45, 0.34, "o", color="k", fontsize=25)

    fig.text(0.64, 0.91, "$\epsilon_\\alpha\equiv-\ln\Lambda_\\alpha^2$", color="k", fontsize=14)

    plt.savefig("/home/bart/Documents/presentations/2021_03_18/figures/phi_flow_nu_1_3_third_new.png", bbox_inches='tight', dpi=300)
    plt.show()