"""Hardcore boson Haldane model.
Hamiltonian based on: "Characterizing topological order by studying the ground states of an infinite cylinder"."""

import numpy as np

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite


class BosonicHaldaneModel(CouplingMPOModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', 'Honeycomb')
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', 1/3, self.name)
        site = BosonSite(conserve=conserve, filling=filling, Nmax=1)
        return site

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        phi_ext = - 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        phi = 0.4*np.pi
        tdash = 0.6*t * np.exp(1j * phi)
        tdashdash = -0.58*t

        for u1, u2, dx in self.lat.nearest_neighbors:

            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext])

            self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
            self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.

        for u1, u2, dx in [(0, 0, np.array([1, 0])), (0, 0, np.array([0, -1])), (0, 0, np.array([-1, 1])),
                           (1, 1, np.array([-1, 0])), (1, 1, np.array([0, 1])), (1, 1, np.array([1, -1]))]:

            tdash_phi = self.coupling_strength_add_ext_flux(tdash, dx, [0, phi_ext])

            self.add_coupling(tdash_phi, u1, 'Bd', u2, 'B', dx)
            self.add_coupling(np.conj(tdash_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.

        for u1, u2, dx in self.lat.next_next_nearest_neighbors:

            tdashdash_phi = self.coupling_strength_add_ext_flux(tdashdash, dx, [0, phi_ext])

            self.add_coupling(tdashdash_phi, u1, 'Bd', u2, 'B', dx)
            self.add_coupling(np.conj(tdashdash_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.


class BosonicHaldaneChain(BosonicHaldaneModel, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
