"""tri_2 model"""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
from lattices.MagneticTriangular import MagneticTriangular


class BosonicTri1Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        Nmax = get_parameter(model_params, 'Nmax', 1, self.name)
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', (1, 8), self.name)
        filling = filling[0] / filling[1]
        site = BosonSite(Nmax=Nmax, conserve=conserve, filling=filling)
        return site

    def init_lattice(self, model_params):
        bc_MPS = get_parameter(model_params, 'bc_MPS', 'infinite', self.name)
        order = get_parameter(model_params, 'order', 'default', self.name)
        site = self.init_sites(model_params)
        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 4, self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = MagneticTriangular(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)

        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 4), self.name)
        phi = 2 * np.pi * phi_p / phi_q

        for i in range(phi_q):
            for u1, u2, dx in getattr(self.lat, "NN_u{}".format(i)):
                # print("lat_NN_u =", u1, u2, dx)
                t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                        * np.exp(1j * (phi/2) * ((2*i)+1/2))
                self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "NN_rd{}".format(i)):
                # print("lat_NN_ul =", u1, u2, dx)
                t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi/2) * ((2*i)+1/2))
                self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "NN_ru{}".format(i)):
                # print("lat_NN_ur =", u1, u2, dx)
                self.add_coupling(t, u1, 'Bd', u2, 'B', dx)
                self.add_coupling(np.conj(t), u2, 'Bd', u1, 'B', -dx)  # h.c.


class FermionicTri1Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', (1, 9), self.name)
        filling = filling[0] / filling[1]
        site = FermionSite(conserve=conserve, filling=filling)
        return site

    def init_lattice(self, model_params):
        bc_MPS = get_parameter(model_params, 'bc_MPS', 'infinite', self.name)
        order = get_parameter(model_params, 'order', 'default', self.name)
        site = self.init_sites(model_params)
        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 6, self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = MagneticTriangular(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        V = get_parameter(model_params, 'V', 10, self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)

        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 3), self.name)
        phi = 2 * np.pi * phi_p / phi_q

        for i in range(phi_q):
            for u1, u2, dx in getattr(self.lat, "NN_u{}".format(i)):
                # print("lat_NN_u =", u1, u2, dx)
                t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                        * np.exp(1j * (phi/2) * ((2*i)+1/2))
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
                self.add_coupling(V, u1, 'N', u2, 'N', dx)
            for u1, u2, dx in getattr(self.lat, "NN_rd{}".format(i)):
                # print("lat_NN_ul =", u1, u2, dx)
                t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi/2) * ((2*i)+1/2))
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
                self.add_coupling(V, u1, 'N', u2, 'N', dx)
            for u1, u2, dx in getattr(self.lat, "NN_ru{}".format(i)):
                # print("lat_NN_ur =", u1, u2, dx)
                self.add_coupling(t, u1, 'Cd', u2, 'C', dx)
                self.add_coupling(np.conj(t), u2, 'Cd', u1, 'C', -dx)  # h.c.
                self.add_coupling(V, u1, 'N', u2, 'N', dx)
