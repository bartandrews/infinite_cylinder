import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# evenly sampled time at 200ms intervals
t = np.arange(0., 5., 0.2)

if __name__ == '__main__':

    totvectors = 13

    fig = plt.figure(figsize=(13.75, 13.75))
    gs = gridspec.GridSpec(3, 2, hspace=0.3)

    ax1 = plt.subplot(gs[0])  # ########################################################################################

    shells = [[] for i in range(totvectors)]
    xcomp = [[] for i in range(totvectors)]
    ycomp = [[] for i in range(totvectors)]
    zcomp = [[] for i in range(totvectors)]
    trace = [[] for i in range(totvectors)]

    with open('CubTable.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            if row[0] + row[1] + row[2] == "{{0 0 1}":
                idx = 0
            elif row[0] + row[1] + row[2] == "{{0 1 -1}":
                idx = 1
            elif row[0] + row[1] + row[2] == "{{0 1 0}":
                idx = 2
            elif row[0] + row[1] + row[2] == "{{0 1 1}":
                idx = 3
            elif row[0] + row[1] + row[2] == "{{1 -1 -1}":
                idx = 4
            elif row[0] + row[1] + row[2] == "{{1 -1 0}":
                idx = 5
            elif row[0] + row[1] + row[2] == "{{1 -1 1}":
                idx = 6
            elif row[0] + row[1] + row[2] == "{{1 0 -1}":
                idx = 7
            elif row[0] + row[1] + row[2] == "{{1 0 0}":
                idx = 8
            elif row[0] + row[1] + row[2] == "{{1 0 1}":
                idx = 9
            elif row[0] + row[1] + row[2] == "{{1 1 -1}":
                idx = 10
            elif row[0] + row[1] + row[2] == "{{1 1 0}":
                idx = 11
            elif row[0] + row[1] + row[2] == "{{1 1 1}":
                idx = 12
            shells[idx].append(float(row[3]))
            xcomp[idx].append(float(row[4]))
            ycomp[idx].append(float(row[5]))
            zcomp[idx].append(float(row[6]))
            trace[idx].append(float(row[7].replace("}", "").replace("*^", "E").replace(" ", "")))

    vectors = ["$(0, 0, 1)$", "$(0, 1, -1)$", "$(0, 1, 0)$", "$(0, 1, 1)$", "$(1, -1, -1)$", "$(1, -1, 0)$",
               "$(1, -1, 1)$", "$(1, 0, -1)$", "$(1, 0, 0)$", "$(1, 0, 1)$", "$(1, 1, -1)$", "$(1, 1, 0)$",
               "$(1, 1, 1)$"]
    colors=["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "brown", "darkgreen", "black"]
    for i in range(totvectors):
        ax1.plot(shells[i], xcomp[i], c=colors[i], marker='^', linestyle='--')
        ax1.plot(shells[i], ycomp[i], c=colors[i], marker='s', linestyle='--')
        ax1.plot(shells[i], zcomp[i], c=colors[i], marker='x', linestyle='--')
        ax1.plot(shells[i], trace[i], c=colors[i], label=vectors[i])

    ax1.legend(bbox_to_anchor=(1.1, 1.5), loc="upper center", handletextpad=0.2, columnspacing=0.25, borderpad=0.4,
               framealpha=1, edgecolor='k',
               markerscale=1, fontsize=12, ncol=7, title='$(i_\mathrm{x},i_\mathrm{y},i_\mathrm{z})$', title_fontsize=12)

    ax1.set_xlabel("$n$", fontsize=12)
    ax1.set_ylabel("$\sum_\\alpha \Phi_{Ii\\alpha,Ii\\alpha}$", fontsize=12)
    ax1.tick_params(axis="x", labelsize=11)
    ax1.tick_params(axis="y", labelsize=11)
    # ax1.text(7, 2, "cub", fontsize=14)

    ax2 = plt.subplot(gs[1])  # ########################################################################################

    shells = [[] for i in range(totvectors)]
    xcomp = [[] for i in range(totvectors)]
    ycomp = [[] for i in range(totvectors)]
    zcomp = [[] for i in range(totvectors)]
    trace = [[] for i in range(totvectors)]

    with open('BccTable.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            if row[0] + row[1] + row[2] == "{{0 0 1}":
                idx = 0
            elif row[0] + row[1] + row[2] == "{{0 1 -1}":
                idx = 1
            elif row[0] + row[1] + row[2] == "{{0 1 0}":
                idx = 2
            elif row[0] + row[1] + row[2] == "{{0 1 1}":
                idx = 3
            elif row[0] + row[1] + row[2] == "{{1 -1 -1}":
                idx = 4
            elif row[0] + row[1] + row[2] == "{{1 -1 0}":
                idx = 5
            elif row[0] + row[1] + row[2] == "{{1 -1 1}":
                idx = 6
            elif row[0] + row[1] + row[2] == "{{1 0 -1}":
                idx = 7
            elif row[0] + row[1] + row[2] == "{{1 0 0}":
                idx = 8
            elif row[0] + row[1] + row[2] == "{{1 0 1}":
                idx = 9
            elif row[0] + row[1] + row[2] == "{{1 1 -1}":
                idx = 10
            elif row[0] + row[1] + row[2] == "{{1 1 0}":
                idx = 11
            elif row[0] + row[1] + row[2] == "{{1 1 1}":
                idx = 12
            shells[idx].append(float(row[3]))
            xcomp[idx].append(float(row[4]))
            ycomp[idx].append(float(row[5]))
            zcomp[idx].append(float(row[6]))
            trace[idx].append(float(row[7].replace("}", "").replace("*^", "E").replace(" ", "")))

    vectors = ["$(0, 0, 1)$", "$(0, 1, -1)$", "$(0, 1, 0)$", "$(0, 1, 1)$", "$(1, -1, -1)$", "$(1, -1, 0)$",
               "$(1, -1, 1)$", "$(1, 0, -1)$", "$(1, 0, 0)$", "$(1, 0, 1)$", "$(1, 1, -1)$", "$(1, 1, 0)$",
               "$(1, 1, 1)$"]
    colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "brown", "darkgreen", "black"]
    for i in range(totvectors):
        ax2.plot(shells[i], xcomp[i], c=colors[i], marker='^', linestyle='--')
        ax2.plot(shells[i], ycomp[i], c=colors[i], marker='s', linestyle='--')
        ax2.plot(shells[i], zcomp[i], c=colors[i], marker='x', linestyle='--')
        ax2.plot(shells[i], trace[i], c=colors[i], label=vectors[i])

    ax2.set_xlabel("$n$", fontsize=12)
    ax2.set_ylabel("$\sum_\\alpha \Phi_{Ii\\alpha,Ii\\alpha}$", fontsize=12)
    ax2.tick_params(axis="x", labelsize=11)
    ax2.tick_params(axis="y", labelsize=11)
    # ax2.text(7, 3, "bcc", fontsize=14)

    ax3 = plt.subplot(gs[2])  # ########################################################################################

    shells = [[] for i in range(totvectors)]
    xcomp = [[] for i in range(totvectors)]
    ycomp = [[] for i in range(totvectors)]
    zcomp = [[] for i in range(totvectors)]
    trace = [[] for i in range(totvectors)]

    with open('FccTable.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            if row[0] + row[1] + row[2] == "{{0 0 1}":
                idx = 0
            elif row[0] + row[1] + row[2] == "{{0 1 -1}":
                idx = 1
            elif row[0] + row[1] + row[2] == "{{0 1 0}":
                idx = 2
            elif row[0] + row[1] + row[2] == "{{0 1 1}":
                idx = 3
            elif row[0] + row[1] + row[2] == "{{1 -1 -1}":
                idx = 4
            elif row[0] + row[1] + row[2] == "{{1 -1 0}":
                idx = 5
            elif row[0] + row[1] + row[2] == "{{1 -1 1}":
                idx = 6
            elif row[0] + row[1] + row[2] == "{{1 0 -1}":
                idx = 7
            elif row[0] + row[1] + row[2] == "{{1 0 0}":
                idx = 8
            elif row[0] + row[1] + row[2] == "{{1 0 1}":
                idx = 9
            elif row[0] + row[1] + row[2] == "{{1 1 -1}":
                idx = 10
            elif row[0] + row[1] + row[2] == "{{1 1 0}":
                idx = 11
            elif row[0] + row[1] + row[2] == "{{1 1 1}":
                idx = 12
            shells[idx].append(float(row[3]))
            xcomp[idx].append(float(row[4]))
            ycomp[idx].append(float(row[5]))
            zcomp[idx].append(float(row[6]))
            trace[idx].append(float(row[7].replace("}", "").replace("*^", "E").replace(" ", "")))

    vectors = ["$(0, 0, 1)$", "$(0, 1, -1)$", "$(0, 1, 0)$", "$(0, 1, 1)$", "$(1, -1, -1)$", "$(1, -1, 0)$",
               "$(1, -1, 1)$", "$(1, 0, -1)$", "$(1, 0, 0)$", "$(1, 0, 1)$", "$(1, 1, -1)$", "$(1, 1, 0)$",
               "$(1, 1, 1)$"]
    colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "brown", "darkgreen", "black"]
    for i in range(totvectors):
        ax3.plot(shells[i], xcomp[i], c=colors[i], marker='^', linestyle='--')
        ax3.plot(shells[i], ycomp[i], c=colors[i], marker='s', linestyle='--')
        ax3.plot(shells[i], zcomp[i], c=colors[i], marker='x', linestyle='--')
        ax3.plot(shells[i], trace[i], c=colors[i], label=vectors[i])

    ax3.set_xlabel("$n$", fontsize=12)
    ax3.set_ylabel("$\sum_\\alpha \Phi_{Ii\\alpha,Ii\\alpha}$", fontsize=12)
    ax3.tick_params(axis="x", labelsize=11)
    ax3.tick_params(axis="y", labelsize=11)
    # ax3.text(7, 15, "fcc", fontsize=14)

    ax4 = plt.subplot(gs[3])  # ########################################################################################

    shells = [[] for i in range(totvectors)]
    xcomp = [[] for i in range(totvectors)]
    ycomp = [[] for i in range(totvectors)]
    zcomp = [[] for i in range(totvectors)]
    trace = [[] for i in range(totvectors)]

    with open('DiaTable.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            if row[0] + row[1] + row[2] == "{{0 0 1}":
                idx = 0
            elif row[0] + row[1] + row[2] == "{{0 1 -1}":
                idx = 1
            elif row[0] + row[1] + row[2] == "{{0 1 0}":
                idx = 2
            elif row[0] + row[1] + row[2] == "{{0 1 1}":
                idx = 3
            elif row[0] + row[1] + row[2] == "{{1 -1 -1}":
                idx = 4
            elif row[0] + row[1] + row[2] == "{{1 -1 0}":
                idx = 5
            elif row[0] + row[1] + row[2] == "{{1 -1 1}":
                idx = 6
            elif row[0] + row[1] + row[2] == "{{1 0 -1}":
                idx = 7
            elif row[0] + row[1] + row[2] == "{{1 0 0}":
                idx = 8
            elif row[0] + row[1] + row[2] == "{{1 0 1}":
                idx = 9
            elif row[0] + row[1] + row[2] == "{{1 1 -1}":
                idx = 10
            elif row[0] + row[1] + row[2] == "{{1 1 0}":
                idx = 11
            elif row[0] + row[1] + row[2] == "{{1 1 1}":
                idx = 12
            shells[idx].append(float(row[3]))
            xcomp[idx].append(float(row[4]))
            ycomp[idx].append(float(row[5]))
            zcomp[idx].append(float(row[6]))
            trace[idx].append(float(row[7].replace("}", "").replace("*^", "E").replace(" ", "")))

    vectors = ["$(0, 0, 1)$", "$(0, 1, -1)$", "$(0, 1, 0)$", "$(0, 1, 1)$", "$(1, -1, -1)$", "$(1, -1, 0)$",
               "$(1, -1, 1)$", "$(1, 0, -1)$", "$(1, 0, 0)$", "$(1, 0, 1)$", "$(1, 1, -1)$", "$(1, 1, 0)$",
               "$(1, 1, 1)$"]
    colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "brown", "darkgreen", "black"]
    for i in range(totvectors):
        ax4.plot(shells[i], xcomp[i], c=colors[i], marker='^', linestyle='--')
        ax4.plot(shells[i], ycomp[i], c=colors[i], marker='s', linestyle='--')
        ax4.plot(shells[i], zcomp[i], c=colors[i], marker='x', linestyle='--')
        ax4.plot(shells[i], trace[i], c=colors[i], label=vectors[i])

    ax4.set_xlabel("$n$", fontsize=12)
    ax4.set_ylabel("$\sum_\\alpha \Phi_{Ii\\alpha,Ii\\alpha}$", fontsize=12)
    ax4.tick_params(axis="x", labelsize=11)
    ax4.tick_params(axis="y", labelsize=11)
    # ax4.text(7, 75, "dia", fontsize=14)

    ax5 = plt.subplot(gs[4])  # ########################################################################################

    shells = [[] for i in range(totvectors)]
    xcomp = [[] for i in range(totvectors)]
    ycomp = [[] for i in range(totvectors)]
    zcomp = [[] for i in range(totvectors)]
    trace = [[] for i in range(totvectors)]

    with open('HcpTable.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            if row[0] + row[1] + row[2] == "{{0 0 1}":
                idx = 0
            elif row[0] + row[1] + row[2] == "{{0 1 -1}":
                idx = 1
            elif row[0] + row[1] + row[2] == "{{0 1 0}":
                idx = 2
            elif row[0] + row[1] + row[2] == "{{0 1 1}":
                idx = 3
            elif row[0] + row[1] + row[2] == "{{1 -1 -1}":
                idx = 4
            elif row[0] + row[1] + row[2] == "{{1 -1 0}":
                idx = 5
            elif row[0] + row[1] + row[2] == "{{1 -1 1}":
                idx = 6
            elif row[0] + row[1] + row[2] == "{{1 0 -1}":
                idx = 7
            elif row[0] + row[1] + row[2] == "{{1 0 0}":
                idx = 8
            elif row[0] + row[1] + row[2] == "{{1 0 1}":
                idx = 9
            elif row[0] + row[1] + row[2] == "{{1 1 -1}":
                idx = 10
            elif row[0] + row[1] + row[2] == "{{1 1 0}":
                idx = 11
            elif row[0] + row[1] + row[2] == "{{1 1 1}":
                idx = 12
            shells[idx].append(float(row[3]))
            xcomp[idx].append(float(row[4]))
            ycomp[idx].append(float(row[5]))
            zcomp[idx].append(float(row[6]))
            trace[idx].append(float(row[7].replace("}", "").replace("*^", "E").replace(" ", "")))

    vectors = ["$(0, 0, 1)$", "$(0, 1, -1)$", "$(0, 1, 0)$", "$(0, 1, 1)$", "$(1, -1, -1)$", "$(1, -1, 0)$",
               "$(1, -1, 1)$", "$(1, 0, -1)$", "$(1, 0, 0)$", "$(1, 0, 1)$", "$(1, 1, -1)$", "$(1, 1, 0)$",
               "$(1, 1, 1)$"]
    colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "brown", "darkgreen", "black"]
    for i in range(totvectors):
        ax5.plot(shells[i], xcomp[i], c=colors[i], marker='^', linestyle='--')
        ax5.plot(shells[i], ycomp[i], c=colors[i], marker='s', linestyle='--')
        ax5.plot(shells[i], zcomp[i], c=colors[i], marker='x', linestyle='--')
        ax5.plot(shells[i], trace[i], c=colors[i], label=vectors[i])

    ax5.set_xlabel("$n$", fontsize=12)
    ax5.set_ylabel("$\sum_\\alpha \Phi_{Ii\\alpha,Ii\\alpha}$", fontsize=12)
    ax5.tick_params(axis="x", labelsize=11)
    ax5.tick_params(axis="y", labelsize=11)
    # ax5.text(7, 8, "hcp", fontsize=14)

    ax6 = plt.subplot(gs[5])  # ########################################################################################

    shells = [[] for i in range(totvectors)]
    xcomp = [[] for i in range(totvectors)]
    ycomp = [[] for i in range(totvectors)]
    zcomp = [[] for i in range(totvectors)]
    trace = [[] for i in range(totvectors)]

    with open('DhcpTable.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            if row[0] + row[1] + row[2] == "{{0 0 1}":
                idx = 0
            elif row[0] + row[1] + row[2] == "{{0 1 -1}":
                idx = 1
            elif row[0] + row[1] + row[2] == "{{0 1 0}":
                idx = 2
            elif row[0] + row[1] + row[2] == "{{0 1 1}":
                idx = 3
            elif row[0] + row[1] + row[2] == "{{1 -1 -1}":
                idx = 4
            elif row[0] + row[1] + row[2] == "{{1 -1 0}":
                idx = 5
            elif row[0] + row[1] + row[2] == "{{1 -1 1}":
                idx = 6
            elif row[0] + row[1] + row[2] == "{{1 0 -1}":
                idx = 7
            elif row[0] + row[1] + row[2] == "{{1 0 0}":
                idx = 8
            elif row[0] + row[1] + row[2] == "{{1 0 1}":
                idx = 9
            elif row[0] + row[1] + row[2] == "{{1 1 -1}":
                idx = 10
            elif row[0] + row[1] + row[2] == "{{1 1 0}":
                idx = 11
            elif row[0] + row[1] + row[2] == "{{1 1 1}":
                idx = 12
            shells[idx].append(float(row[3]))
            xcomp[idx].append(float(row[4]))
            ycomp[idx].append(float(row[5]))
            zcomp[idx].append(float(row[6]))
            trace[idx].append(float(row[7].replace("}", "").replace("*^", "E").replace(" ", "")))

    vectors = ["$(0, 0, 1)$", "$(0, 1, -1)$", "$(0, 1, 0)$", "$(0, 1, 1)$", "$(1, -1, -1)$", "$(1, -1, 0)$",
               "$(1, -1, 1)$", "$(1, 0, -1)$", "$(1, 0, 0)$", "$(1, 0, 1)$", "$(1, 1, -1)$", "$(1, 1, 0)$",
               "$(1, 1, 1)$"]
    colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "brown", "darkgreen", "black"]
    for i in range(totvectors):
        ax6.plot(shells[i], xcomp[i], c=colors[i], marker='^', linestyle='--')
        ax6.plot(shells[i], ycomp[i], c=colors[i], marker='s', linestyle='--')
        ax6.plot(shells[i], zcomp[i], c=colors[i], marker='x', linestyle='--')
        ax6.plot(shells[i], trace[i], c=colors[i], label=vectors[i])

    ax6.set_xlabel("$n$", fontsize=12)
    ax6.set_ylabel("$\sum_\\alpha \Phi_{Ii\\alpha,Ii\\alpha}$", fontsize=12)
    ax6.tick_params(axis="x", labelsize=11)
    ax6.tick_params(axis="y", labelsize=11)
    # ax6.text(7, 30, "dhcp", fontsize=14)

    fig.text(0.08, 0.88, "$\mathrm{(a)}$", fontsize=14)
    fig.text(0.08, 0.6, "$\mathrm{(c)}$", fontsize=14)
    fig.text(0.08, 0.32, "$\mathrm{(e)}$", fontsize=14)

    fig.text(0.5, 0.88, "$\mathrm{(b)}$", fontsize=14)
    fig.text(0.5, 0.6, "$\mathrm{(d)}$", fontsize=14)
    fig.text(0.5, 0.32, "$\mathrm{(f)}$", fontsize=14)

    fig.text(0.427, 0.825, "cub", fontsize=14)
    fig.text(0.85, 0.825, "bcc", fontsize=14)
    fig.text(0.427, 0.54, "fcc", fontsize=14)
    fig.text(0.85, 0.54, "dia", fontsize=14)
    fig.text(0.427, 0.265, "hcp", fontsize=14)
    fig.text(0.85, 0.265, "dhcp", fontsize=14)

    plt.savefig("/home/bart/Documents/papers/QLI001/figures/trace.png", bbox_inches='tight', dpi=300)
    plt.show()