"""Cold atomic (Harper-)Hofstadter model on a strip or cylinder.

"""
# Copyright 2018 TeNPy Developers

import numpy as np
import warnings
import sys

from tenpy.models.lattice import Square
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.model import CouplingModel, MPOModel, CouplingMPOModel
from tenpy.tools.params import get_parameter, unused_parameters
from lattices.MagneticSquare import MagneticSquare

__all__ = ['HofstadterBosons', 'gauge_hopping']


def gauge_hopping(model_params):

    t = get_parameter(model_params, 't', 1., 'Gauge hopping')
    phi_p, phi_q = get_parameter(model_params, 'phi', (1, 3), 'Gauge hopping')
    phi = 2 * np.pi * phi_p / phi_q

    hop_x = -t
    hop_y = -t * np.exp(1.j * phi * np.arange(phi_q)[:, np.newaxis])  # has shape (mx, 1)

    print(np.arange(phi_q))
    print(np.arange(phi_q)[:, np.newaxis])
    print("hop_y = ", hop_y)
    print(hop_y[0][0])

    return hop_x, hop_y


class HofstadterBosons(CouplingModel, MPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        Nmax = get_parameter(model_params, 'Nmax', 1, self.__class__)
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', (1, 8), self.name)
        filling = filling[0] / filling[1]
        site = BosonSite(Nmax=Nmax, conserve=conserve, filling=filling)
        return site

    def init_lattice(self, model_params):

        bc_MPS = get_parameter(model_params, 'bc_MPS', 'infinite', self.name)
        order = get_parameter(model_params, 'order', 'default', self.name)
        site = self.init_sites(model_params)
        Lx = get_parameter(model_params, 'Lx', 4, self.name)
        Ly = get_parameter(model_params, 'Ly', 4, self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = MagneticSquare(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
        return lat

    def init_terms(self, model_params):

        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)
        hop_x, hop_y = gauge_hopping(model_params)

        # dx = np.array([1, 0])
        #
        # self.add_coupling(hop_x, 0, 'Bd', 0, 'B', dx)
        # self.add_coupling(np.conj(hop_x), 0, 'Bd', 0, 'B', -dx)  # h.c.
        #
        # dy = np.array([0, 1])
        #
        # print("hop_y = ", hop_y)
        #
        # hop_y = self.coupling_strength_add_ext_flux(hop_y, dy, [0, phi_ext])
        # self.add_coupling(hop_y, 0, 'Bd', 0, 'B', dy)
        # self.add_coupling(np.conj(hop_y), 0, 'Bd', 0, 'B', -dy)  # h.c.

        t = get_parameter(model_params, 't', 1., 'Gauge hopping')
        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 3), 'Gauge hopping')
        phi = 2 * np.pi * phi_p / phi_q

        for u1, u2, dx in self.lat.NN_horiz:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext])
            self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
            self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.

        for i in range(4):
            for u1, u2, dx in getattr(self.lat, "NN_v{}".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j * phi * i)
                self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.
