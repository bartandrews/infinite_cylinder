"""Honeycomb lattice in a perpendicular magnetic field"""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site
import matplotlib.pyplot as plt

import functions as f

class MagneticHoneycomb(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        numb_sites = int(2*f.qval)

        basis = np.array(([(numb_sites/4) * np.sqrt(3), numb_sites/4], [0., 1]))
        delta = np.array([1 / (2. * np.sqrt(3.)), 0.5])

        pos_list = []
        for i in range(0, int(numb_sites/2)):
            pos_list.append([(i/2) * np.sqrt(3), i/2] - delta / 2.)  # A site
            pos_list.append([(i/2) * np.sqrt(3), i/2] + delta / 2.)  # B site
        pos = tuple(pos_list)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA] * numb_sites, **kwargs)

        # NN ###
        for i in range(0, numb_sites, 2):

            setattr(self, "NN{}d".format(i), [(i, i+1, np.array([0, -1]))])

            if i == 0:
                setattr(self, "NN0ul", [(0, numb_sites-1, np.array([-1, 0]))])
            else:
                setattr(self, "NN{}ul".format(i), [(i, i-1, np.array([0, 0]))])

            setattr(self, "NN{}ur".format(i), [(i, i+1, np.array([0, 0]))])

        # secondNN ###
        for i in range(0, numb_sites):

            setattr(self, "secondNN{}bl".format(i), [(i, i, np.array([0, -1]))])

            if i >= (numb_sites-2):
                setattr(self, "secondNN{}r".format(i), [(i, i-(numb_sites-2), np.array([1, 0]))])
            else:
                setattr(self, "secondNN{}r".format(i), [(i, i+2, np.array([0, 0]))])

            if i < 2:
                setattr(self, "secondNN{}ul".format(i), [(i, i+(numb_sites-2), np.array([-1, 1]))])
            else:
                setattr(self, "secondNN{}ul".format(i), [(i, i-2, np.array([0, 1]))])

        # thirdNN ###
        for i in range(0, numb_sites, 2):

            if i == 0:
                setattr(self, "thirdNN{}u".format(i), [(i, i+(numb_sites-1), np.array([-1, 1]))])
            else:
                setattr(self, "thirdNN{}u".format(i), [(i, i-1, np.array([0, 1]))])

            if i >= (numb_sites-3):
                setattr(self, "thirdNN{}br".format(i), [(i, i - (numb_sites-3), np.array([1, -1]))])
            else:
                setattr(self, "thirdNN{}br".format(i), [(i, i + 3, np.array([0, -1]))])

            if i == 0:
                setattr(self, "thirdNN{}bl".format(i), [(i, i + (numb_sites-1), np.array([-1, -1]))])
            else:
                setattr(self, "thirdNN{}bl".format(i), [(i, i -1, np.array([0, -1]))])

        # fourthNN ###
        for i in range(0, numb_sites, 2):

            if i >= (numb_sites-3):
                setattr(self, "fourthNN{}bur".format(i), [(i, i - (numb_sites-3), np.array([1, 0]))])
            else:
                setattr(self, "fourthNN{}bur".format(i), [(i, i + 3, np.array([0, 0]))])

            if i >= (numb_sites-1):
                setattr(self, "fourthNN{}uur".format(i), [(i, i - (numb_sites-1), np.array([1, 1]))])
            else:
                setattr(self, "fourthNN{}uur".format(i), [(i, i + 1, np.array([0, 1]))])

            if i < 3:
                setattr(self, "fourthNN{}bul".format(i), [(i, i + (numb_sites-3), np.array([-1, 0]))])
            else:
                setattr(self, "fourthNN{}bul".format(i), [(i, i - 3, np.array([0, 0]))])

            if i < 3:
                setattr(self, "fourthNN{}uul".format(i), [(i, i + (numb_sites-3), np.array([-1, 1]))])
            else:
                setattr(self, "fourthNN{}uul".format(i), [(i, i - 3, np.array([0, 1]))])

            if i >= (numb_sites-1):
                setattr(self, "fourthNN{}bl".format(i), [(i, i - (numb_sites-1), np.array([1, -2]))])
            else:
                setattr(self, "fourthNN{}bl".format(i), [(i, i + 1, np.array([0, -2]))])

            if i >= (numb_sites-3):
                setattr(self, "fourthNN{}br".format(i), [(i, i - (numb_sites-3), np.array([1, -2]))])
            else:
                setattr(self, "fourthNN{}br".format(i), [(i, i + 3, np.array([0, -2]))])

        # fifthNN ###
        for i in range(0, numb_sites):

            if i < 2:
                setattr(self, "fifthNN{}u".format(i), [(i, i+(numb_sites-2), np.array([-1, 2]))])
            else:
                setattr(self, "fifthNN{}u".format(i), [(i, i-2, np.array([0, 2]))])

            if i < 2:
                setattr(self, "fifthNN{}bl".format(i), [(i, i+(numb_sites-2), np.array([-1, -1]))])
            else:
                setattr(self, "fifthNN{}bl".format(i), [(i, i-2, np.array([0, -1]))])

            if i >= (numb_sites-4):
                setattr(self, "fifthNN{}br".format(i), [(i, i-(numb_sites-4), np.array([1, -1]))])
            else:
                setattr(self, "fifthNN{}br".format(i), [(i, i+4, np.array([0, -1]))])


def plot_lattice():

    numb_sites = int(2*f.qval)

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticHoneycomb(1, 4, fs, basis=[[(numb_sites/4) * np.sqrt(3), numb_sites/4], [0, 1]])
    lat.plot_sites(ax)

    # for i in range(0, numb_sites, 2):
    #     lat.plot_coupling(ax, getattr(lat, "NN{}d".format(i)), linestyle='-', color='red')
    #     lat.plot_coupling(ax, getattr(lat, "NN{}ul".format(i)), linestyle='-', color='blue')
    #     lat.plot_coupling(ax, getattr(lat, "NN{}ur".format(i)), linestyle='-', color='green')
    for i in range(0, numb_sites):
        lat.plot_coupling(ax, getattr(lat, "secondNN{}bl".format(i)), linestyle='-', color='red')
        lat.plot_coupling(ax, getattr(lat, "secondNN{}r".format(i)), linestyle='-', color='blue')
        lat.plot_coupling(ax, getattr(lat, "secondNN{}ul".format(i)), linestyle='-', color='green')
    # for i in range(0, numb_sites, 2):
    #     lat.plot_coupling(ax, getattr(lat, "thirdNN{}u".format(i)), linestyle='-', color='red')
    #     lat.plot_coupling(ax, getattr(lat, "thirdNN{}br".format(i)), linestyle='-', color='blue')
    #     lat.plot_coupling(ax, getattr(lat, "thirdNN{}bl".format(i)), linestyle='-', color='green')
    # for i in range(0, numb_sites, 2):
    #     lat.plot_coupling(ax, getattr(lat, "fourthNN{}bur".format(i)), linestyle='-', color='red')
    #     lat.plot_coupling(ax, getattr(lat, "fourthNN{}uur".format(i)), linestyle='-', color='blue')
    #     lat.plot_coupling(ax, getattr(lat, "fourthNN{}bul".format(i)), linestyle='-', color='green')
    #     lat.plot_coupling(ax, getattr(lat, "fourthNN{}uul".format(i)), linestyle='-', color='orange')
    #     lat.plot_coupling(ax, getattr(lat, "fourthNN{}bl".format(i)), linestyle='-', color='purple')
    #     lat.plot_coupling(ax, getattr(lat, "fourthNN{}br".format(i)), linestyle='-', color='black')
    # for i in range(0, numb_sites):
    #     lat.plot_coupling(ax, getattr(lat, "fifthNN{}u".format(i)), linestyle='-', color='orange')
    #     lat.plot_coupling(ax, getattr(lat, "fifthNN{}bl".format(i)), linestyle='-', color='purple')
    #     lat.plot_coupling(ax, getattr(lat, "fifthNN{}br".format(i)), linestyle='-', color='black')

    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":

    plot_lattice()
