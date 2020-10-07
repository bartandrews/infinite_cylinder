import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

if __name__ == '__main__':

    fig = plt.figure(figsize=(6, 9))
    gs = gridspec.GridSpec(5, 1, hspace=0.7)

    # bifurcation plot #################################################################################################

    # initialization
    dr, dx = 0.01, 0.01
    max_it = 10000
    R = np.arange(0, 4, dr)
    X = np.arange(0, 1, dx)
    Z = np.zeros((len(X), len(R)))
    lyapunov = np.zeros((max_it, len(R)))  # for lyapunov exponent plot
    zarray = np.zeros((max_it, len(R)))  # for logistic plot

    # evaluation
    for ir, r in enumerate(R):
        if ir % 100 == 0:
            print(ir, "of", len(R))
        z = 0.25  # any value in [0,1]
        for i in range(1000):  # stabilization (throw away initial behaviour)
            z = r * z * (1 - z)
        for i in range(max_it):  # iteration (examine asympototic behaviour)
            z = r * z * (1 - z)
            Z[int(z / dx), ir] += 1
            zarray[i, ir] = z
            lyapunov[i, ir] = np.log(abs(r * (1 - 2 * z)))
    Z = np.where(Z > 0, 1, np.NaN)

    # plotting
    ax = plt.subplot(gs[0])
    ax.imshow(Z, interpolation='none', origin='lower', extent=(min(R), max(R), min(X), max(X)), aspect='auto')
    ax.set_xlabel("$r$")
    ax.set_xlim([0, 4])
    ax.set_xticks(np.arange(0, 4.1, 0.5))
    ax.set_ylabel("$x_{n+1}=f(x_n)$")
    ax.set_ylim([0, 1])
    ax.set_yticks(np.arange(0, 1.1, 0.5))

    # lyapunov exponent plot ###########################################################################################

    ax1 = plt.subplot(gs[1], sharex=ax)
    ax1.plot(R, lyapunov.mean(axis=0))
    ax1.set_xlabel("$r$")
    ax1.set_xlim([0, 4])
    ax1.set_xticks(np.arange(0, 4.1, 0.5))
    ax1.set_ylabel("$\lambda$")
    ax1.axhline(0, color='k', linewidth=0.5, ls='--')

    # shannon entropy ##################################################################################################

    decimal_places = 2

    shannonarray = np.zeros((max_it, len(R)))
    shannon = np.zeros(len(R))

    for ir, r in enumerate(R):
        unique, counts = np.unique(np.round(zarray[:, ir], decimals=decimal_places), return_counts=True)
        dic = dict(zip(unique, counts))
        for i in range(len(dic)):
            p = dic[np.round(zarray[i, ir], decimals=decimal_places)] / max_it
            shannonarray[i, ir] = - p * np.log2(p)
        shannon[ir] = sum(shannonarray[:, ir])

    ax2 = plt.subplot(gs[2], sharex=ax)
    ax2.plot(R, shannon, 'x', markersize=1)
    ax2.set_xlabel("$r$")
    ax2.set_xlim([0, 4])
    ax2.set_xticks(np.arange(0, 4.1, 0.5))
    ax2.set_ylabel("$S$")

    # kolmogorov entropy ###############################################################################################

    # positive Lyapunov exponents

    kolmogorov = lyapunov.mean(axis=0)
    kolmogorov[kolmogorov < 0] = 0

    # decimal_places = 2
    # N = 10
    #
    # difference in shannon entropies of different length orbits
    #
    # shannon_collec = np.zeros((N, len(R)))
    # kolmogorovarray = np.zeros((N - 1, len(R)))
    # kolmogorov = np.zeros(len(R))
    #
    # for j in range(1, N):
    #
    #     shannonarray = np.zeros((j, len(R)))
    #     shannon = np.zeros(len(R))
    #
    #     for ir, r in enumerate(R):
    #         unique, counts = np.unique(np.round(zarray[:j, ir], decimals=decimal_places), return_counts=True)
    #         dic = dict(zip(unique, counts))
    #         for i in range(len(dic)):
    #             p = dic[np.round(zarray[i, ir], decimals=decimal_places)] / j
    #             shannonarray[i, ir] = - p * np.log2(p)
    #         shannon[ir] = sum(shannonarray[:, ir])
    #         if ir == 305:
    #             print(ir, j, shannon[ir])
    #
    #         shannon_collec[j-1, ir] = shannon[ir]
    #         if j > 1:
    #             kolmogorovarray[j-1, ir] = shannon_collec[j-1, ir] - shannon_collec[j - 2, ir]
    #             if ir == 305:
    #                 print(kolmogorovarray[j-1, ir])
    #
    # for ir, r in enumerate(R):
    #     kolmogorov[ir] = np.mean(kolmogorovarray[:, ir])

    # difference in shannon entropies of same length orbits
    #
    # kolmogorovarray = np.zeros((max_it-1, len(R)))
    # kolmogorov = np.zeros(len(R))
    # for ir, r in enumerate(R):
    #     unique, counts = np.unique(np.round(zarray[:, ir], decimals=decimal_places), return_counts=True)
    #     dic = dict(zip(unique, counts))
    #     for i in range(len(dic)-1):
    #         p_n = dic[np.round(zarray[i, ir], decimals=decimal_places)] / max_it
    #         p_n_plus_1 = dic[np.round(zarray[i+1, ir], decimals=decimal_places)] / max_it
    #         kolmogorovarray[i, ir] = - p_n_plus_1 * np.log2(p_n_plus_1) + p_n * np.log2(p_n)
    #     kolmogorov[ir] = abs(sum(kolmogorovarray[:, ir]))

    ax3 = plt.subplot(gs[3], sharex=ax)
    ax3.plot(R, kolmogorov, 'x', markersize=1)
    ax3.set_xlabel("$r$")
    ax3.set_xlim([0, 4])
    ax3.set_xticks(np.arange(0, 4.1, 0.5))
    ax3.set_ylabel("$K$")

    # logistic plot ####################################################################################################

    rvalue = 3.2

    ax4 = plt.subplot(gs[4])
    ax4.plot(zarray[:100, int(rvalue*100)], linewidth=0.5, label=f"$r={rvalue}$")
    ax4.legend(loc='upper right')
    ax4.set_xlabel("$n$")
    ax4.set_xlim([0, 100])
    ax4.set_xticks(np.arange(0, 101, 10))
    ax4.set_ylabel("$x_n$")
    ax.axvline(rvalue, color='k', linewidth=0.5, ls='--')

    ####################################################################################################################

    plt.suptitle("$x_{n+1}=r x_n (1 - x_n)$", y=0.95)
    plt.savefig("/home/bart/Documents/papers/SR/figures/logistic_map.png", bbox_inches='tight', dpi=300)
    plt.show()
