"""Twist lattice in a perpendicular magnetic field (Landau gauge in x direction with alpha=2/7)."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MagneticTwist(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        basis = np.array(([3.5 * np.sqrt(3), 3.5], [0., 1]))
        delta = np.array([1 / (2. * np.sqrt(3.)), 0.5])
        pos = (-delta / 2., delta / 2, [0.5 * np.sqrt(3), 0.5] - delta / 2., [0.5 * np.sqrt(3), 0.5] + delta / 2,
               [1 * np.sqrt(3), 1] - delta / 2., [1 * np.sqrt(3), 1] + delta / 2, [1.5 * np.sqrt(3), 1.5] - delta / 2., [1.5 * np.sqrt(3), 1.5] + delta / 2,
               [2 * np.sqrt(3), 2] - delta / 2., [2 * np.sqrt(3), 2] + delta / 2, [2.5 * np.sqrt(3), 2.5] - delta / 2., [2.5 * np.sqrt(3), 2.5] + delta / 2,
               [3 * np.sqrt(3), 3] - delta / 2., [3 * np.sqrt(3), 3] + delta / 2)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA, siteA, siteA, siteA, siteA, siteA,
                                    siteA, siteA, siteA, siteA, siteA, siteA, siteA, siteA], **kwargs)

        self.NN0d = [(0, 1, np.array([0, -1]))]
        self.NN2d = [(2, 3, np.array([0, -1]))]
        self.NN4d = [(4, 5, np.array([0, -1]))]
        self.NN6d = [(6, 7, np.array([0, -1]))]
        self.NN8d = [(8, 9, np.array([0, -1]))]
        self.NN10d = [(10, 11, np.array([0, -1]))]
        self.NN12d = [(12, 13, np.array([0, -1]))]

        self.NN0ul = [(0, 9, np.array([-1, 0]))]
        self.NN2ul = [(2, 1, np.array([0, 0]))]
        self.NN4ul = [(4, 3, np.array([0, 0]))]
        self.NN6ul = [(6, 5, np.array([0, 0]))]
        self.NN8ul = [(8, 7, np.array([0, 0]))]
        self.NN10ul = [(10, 9, np.array([0, 0]))]
        self.NN12ul = [(12, 11, np.array([0, 0]))]

        self.NN0ur = [(0, 1, np.array([0, 0]))]
        self.NN2ur = [(2, 3, np.array([0, 0]))]
        self.NN4ur = [(4, 5, np.array([0, 0]))]
        self.NN6ur = [(6, 7, np.array([0, 0]))]
        self.NN8ur = [(8, 9, np.array([0, 0]))]
        self.NN10ur = [(10, 11, np.array([0, 0]))]
        self.NN12ur = [(12, 13, np.array([0, 0]))]

        self.fifthNN0u = [(0, 12, np.array([-1, 2]))]
        self.fifthNN1u = [(1, 13, np.array([-1, 2]))]
        self.fifthNN2u = [(2, 0, np.array([0, 2]))]
        self.fifthNN3u = [(3, 1, np.array([0, 2]))]
        self.fifthNN4u = [(4, 2, np.array([0, 2]))]
        self.fifthNN5u = [(5, 3, np.array([0, 2]))]
        self.fifthNN6u = [(6, 4, np.array([0, 2]))]
        self.fifthNN7u = [(7, 5, np.array([0, 2]))]
        self.fifthNN8u = [(8, 6, np.array([0, 2]))]
        self.fifthNN9u = [(9, 7, np.array([0, 2]))]
        self.fifthNN10u = [(10, 8, np.array([0, 2]))]
        self.fifthNN11u = [(11, 9, np.array([0, 2]))]
        self.fifthNN12u = [(12, 10, np.array([0, 2]))]
        self.fifthNN13u = [(13, 11, np.array([0, 2]))]

        self.fifthNN0bl = [(0, 12, np.array([-1, -1]))]
        self.fifthNN1bl = [(1, 13, np.array([-1, -1]))]
        self.fifthNN2bl = [(2, 0, np.array([0, -1]))]
        self.fifthNN3bl = [(3, 1, np.array([0, -1]))]
        self.fifthNN4bl = [(4, 2, np.array([0, -1]))]
        self.fifthNN5bl = [(5, 3, np.array([0, -1]))]
        self.fifthNN6bl = [(6, 4, np.array([0, -1]))]
        self.fifthNN7bl = [(7, 5, np.array([0, -1]))]
        self.fifthNN8bl = [(8, 6, np.array([0, -1]))]
        self.fifthNN9bl = [(9, 7, np.array([0, -1]))]
        self.fifthNN10bl = [(10, 8, np.array([0, -1]))]
        self.fifthNN11bl = [(11, 9, np.array([0, -1]))]
        self.fifthNN12bl = [(12, 10, np.array([0, -1]))]
        self.fifthNN13bl = [(13, 11, np.array([0, -1]))]

        self.fifthNN0br = [(0, 4, np.array([0, -1]))]
        self.fifthNN1br = [(1, 5, np.array([0, -1]))]
        self.fifthNN2br = [(2, 6, np.array([0, -1]))]
        self.fifthNN3br = [(3, 7, np.array([0, -1]))]
        self.fifthNN4br = [(4, 8, np.array([0, -1]))]
        self.fifthNN5br = [(5, 9, np.array([0, -1]))]
        self.fifthNN6br = [(6, 10, np.array([0, -1]))]
        self.fifthNN7br = [(7, 11, np.array([0, -1]))]
        self.fifthNN8br = [(8, 12, np.array([0, -1]))]
        self.fifthNN9br = [(9, 13, np.array([0, -1]))]
        self.fifthNN10br = [(10, 0, np.array([1, -1]))]
        self.fifthNN11br = [(11, 1, np.array([1, -1]))]
        self.fifthNN12br = [(12, 2, np.array([1, -1]))]
        self.fifthNN13br = [(13, 3, np.array([1, -1]))]


def plot_lattice():

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticTwist(5, 5, fs, basis=[[3.5 *np.sqrt(3), 3.5], [0, 1]])
    lat.plot_sites(ax)
    # ## NN
    # # down
    # lat.plot_coupling(ax, lat.NN0d, linestyle='-', color='black')
    # lat.plot_coupling(ax, lat.NN2d, linestyle='-', color='red')
    # lat.plot_coupling(ax, lat.NN4d, linestyle='-', color='blue')
    # lat.plot_coupling(ax, lat.NN6d, linestyle='-', color='green')
    # lat.plot_coupling(ax, lat.NN8d, linestyle='-', color='orange')
    # lat.plot_coupling(ax, lat.NN10d, linestyle='-', color='black')
    # lat.plot_coupling(ax, lat.NN12d, linestyle='-', color='red')
    # # upper left
    # lat.plot_coupling(ax, lat.NN0ul, linestyle='-', color='black')
    # lat.plot_coupling(ax, lat.NN2ul, linestyle='-', color='red')
    # lat.plot_coupling(ax, lat.NN4ul, linestyle='-', color='blue')
    # lat.plot_coupling(ax, lat.NN6ul, linestyle='-', color='green')
    # lat.plot_coupling(ax, lat.NN8ul, linestyle='-', color='orange')
    # lat.plot_coupling(ax, lat.NN10ul, linestyle='-', color='black')
    # lat.plot_coupling(ax, lat.NN12ul, linestyle='-', color='red')
    # # upper right
    # lat.plot_coupling(ax, lat.NN0ur, linestyle='-', color='black')
    # lat.plot_coupling(ax, lat.NN2ur, linestyle='-', color='red')
    # lat.plot_coupling(ax, lat.NN4ur, linestyle='-', color='blue')
    # lat.plot_coupling(ax, lat.NN6ur, linestyle='-', color='green')
    # lat.plot_coupling(ax, lat.NN8ur, linestyle='-', color='orange')
    # lat.plot_coupling(ax, lat.NN10ur, linestyle='-', color='black')
    # lat.plot_coupling(ax, lat.NN12ur, linestyle='-', color='red')
    ## fifthNN
    # up
    lat.plot_coupling(ax, lat.fifthNN0u, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN1u, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN2u, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN3u, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.fifthNN4u, linestyle='-', color='orange')
    lat.plot_coupling(ax, lat.fifthNN5u, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN6u, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN7u, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN8u, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.fifthNN9u, linestyle='-', color='orange')
    lat.plot_coupling(ax, lat.fifthNN10u, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN11u, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN12u, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN13u, linestyle='-', color='green')
    # bottom left
    lat.plot_coupling(ax, lat.fifthNN0bl, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN1bl, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN2bl, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN3bl, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.fifthNN4bl, linestyle='-', color='orange')
    lat.plot_coupling(ax, lat.fifthNN5bl, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN6bl, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN7bl, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN8bl, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.fifthNN9bl, linestyle='-', color='orange')
    lat.plot_coupling(ax, lat.fifthNN10bl, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN11bl, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN12bl, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN13bl, linestyle='-', color='green')
    # bottom right
    lat.plot_coupling(ax, lat.fifthNN0br, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN1br, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN2br, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN3br, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.fifthNN4br, linestyle='-', color='orange')
    lat.plot_coupling(ax, lat.fifthNN5br, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN6br, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN7br, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN8br, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.fifthNN9br, linestyle='-', color='orange')
    lat.plot_coupling(ax, lat.fifthNN10br, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.fifthNN11br, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.fifthNN12br, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.fifthNN13br, linestyle='-', color='green')
    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
