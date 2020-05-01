import numpy as np
import matplotlib.pyplot as plt


def cost(q_val, Ly_val, q_min_val, Ly_min_val):
    cost_val = q_val * np.exp(Ly_val) / (q_min_val*np.exp(Ly_min_val))  # normalization factor
    return cost_val


def LylB(nphi_val, Ly_val):
    return np.sqrt(2*np.pi*nphi_val)*Ly_val


if __name__ == '__main__':

    LylB_separation = 0.5  # keep all LylB values at least this distance away from each other

    counter = 0
    for Ly in range(3, 16):  # adjust 3 <= Ly <= 15
        for p in range(1, 21):
            for q in range(4, 21):
                nphi = p/q
                if (10/Ly)**2/(2*np.pi) < nphi < np.minimum(0.4, (15/Ly)**2/(2*np.pi)):  # adjust 10 < LylB < 15
                    if counter == 0:
                        data = np.array([[Ly, p, q, LylB(nphi, Ly), 1]])
                        q_min, Ly_min = q, Ly
                    else:
                        if all(abs(i - LylB(nphi, Ly)) >= LylB_separation for i in list(data[:, 3])):
                            data_line = np.array([[Ly, p, q, LylB(nphi, Ly), cost(q, Ly, q_min, Ly_min)]])
                            data = np.concatenate((data, data_line))
                    counter += 1

    # sort the array by cost
    sorted_array = data[np.argsort(data[:, 4])]

    # print to the screen
    print("Ly\tp\tq\tLylB\tcost")
    for line in sorted_array:
        print("{}\t{}\t{}\t{:.3f}\t{:.3f}".format(int(line[0]), int(line[1]), int(line[2]), LylB(line[1]/line[2], line[0]), line[4]))

    # plot the figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    LylB = [LylB(i[1]/i[2], i[0]) for i in sorted_array]
    y = range(1, len(LylB)+1)
    ax.scatter(LylB, y, c=y)
    ax.set_xlabel("$L_y/l_B$")
    ax.set_ylabel("$n_\phi$ index (sorted by increasing numerical cost)")
    plt.show()
