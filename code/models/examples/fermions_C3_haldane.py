"""Spinless fermion generalized C=3 Haldane model.
Hamiltonian based on: "Topological flat band models with arbitrary Chern numbers"."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite, GroupedSite
from lattices.TripartiteTriangular import TripartiteTriangular


class FermionicC3HaldaneModel(CouplingMPOModel):

    def __init__(self, model_params):

        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        fs = FermionSite(conserve=conserve)

        gs = GroupedSite([fs, fs], labels=['A', 'B'], charges='same')
        gs.add_op('Ntot', gs.NA + gs.NB, False)

        return gs

    def init_lattice(self, model_params):

        choice = get_parameter(model_params, 'lattice', 'TripartiteTriangular', self.name)

        if choice != 'TripartiteTriangular':
            sys.exit("Error: Please choose the TripartiteTriangular for C3_haldane.")

        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 3, self.name)

        fs = self.init_sites(model_params)

        lat = TripartiteTriangular(Lx, Ly, fs)

        print(lat.N_sites)

        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        mu = get_parameter(model_params, 'mu', 0., self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        t1 = t
        t2 = 0.39*t*1j
        t3 = -0.34*t

        for u in range(len(self.lat.unit_cell)):

            self.add_onsite(mu, 0, 'N', category='mu N')
            self.add_onsite(-mu, 1, 'N', category='mu N')

        for u1, u2, dx in self.lat.NN:

            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext])
            self.add_coupling(t1_phi, u1, 'CdA', u2, 'CB', dx, 'JW', True)
            self.add_coupling(np.conj(t1_phi), u2, 'CdB', u1, 'CA', -dx, 'JW', True)

        for u1, u2, dx in self.lat.nNNA:

            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext])
            self.add_coupling(t2_phi, u1, 'CdA', u2, 'CA', dx, 'JW', True)
            self.add_coupling(np.conj(t2_phi), u2, 'CdA', u1, 'CA', -dx, 'JW', True)

        for u1, u2, dx in self.lat.nNNB:

            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext])
            self.add_coupling(t2_phi, u1, 'CdB', u2, 'CB', dx, 'JW', True)
            self.add_coupling(np.conj(t2_phi), u2, 'CdB', u1, 'CB', -dx, 'JW', True)

        for u1, u2, dx in self.lat.nnNN:

            t3_phi = self.coupling_strength_add_ext_flux(t3, dx, [0, phi_ext])
            self.add_coupling(t3_phi, u1, 'CdA', u2, 'CB', dx, 'JW', True)
            self.add_coupling(np.conj(t3_phi), u2, 'CdB', u1, 'CA', -dx, 'JW', True)


class FermionicC3HaldaneChain(FermionicC3HaldaneModel, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
