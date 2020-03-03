"""Hofstadter model.
Hamiltonian based on: "Square Lattice with Magnetic Field", Aidelsburger PhD thesis."""

import numpy as np
import sys
import time

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.lattice import Square


class BosonicHofstadterModel(CouplingMPOModel):
    r"""Hardcore bosonic Haldane model.

        The Hamiltonian reads:

        .. math ::
            H = \sum_{ij} t_{ij} b_i^\dagger b_j + \sum_i \mu (n_{A, i} - n_{B, i})
            + V \sum_{\langle ij \rangle, i<j} n_{A, i} n_{B, j}


        Here, :math:`\langle i,j \rangle, i< j` denotes nearest neighbor pairs and :math:`n_A, n_B` are the number operators
        on the A and B sites. Hopping is allowed to nearest and next-nearest neighbor sites with amplitudes
        :math:`t_{\langle ij \rangle}=t_1 \in \mathbb{R}` and
        :math:`t_{\langle\langle ij \rangle\rangle}=t_2 e^{\pm\mathrm{i}\phi} \in \mathbb{C}` respectively, where
        :math:`\pm\phi` is the phase acquired by a boson hopping between atoms in the same sublattice with a sign
        given by the direction of the hopping. This Hamiltonian is translated from [Grushin2015]_.
        All parameters are collected in a single dictionary `model_params` and read out with
        :func:`~tenpy.tools.params.get_parameter`.

        Parameters
        ----------
        conserve : 'best' | 'N' | 'parity' | None
            What should be conserved. See :class:`~tenpy.networks.Site.BosonSite`.
            For ``'best'``, we check the parameters that can be preserved.
        t1, t2, V, mu : float | array
            Hopping, interaction and chemical potential as defined for the Hamiltonian above.
            The default value for t2 is chosen to achieve the optimal band flatness ratio.
        bc_MPS : {'finite' | 'infinte'}
            MPS boundary conditions along the x-direction.
            For 'infinite' boundary conditions, repeat the unit cell in x-direction.
            Coupling boundary conditions in x-direction are chosen accordingly.
            Only used if `lattice` is a string.
        order : string
            Ordering of the sites in the MPS, e.g. 'default', 'snake';
            see :meth:`~tenpy.models.lattice.Lattice.ordering`.
            Only used if `lattice` is a string.
        L : int
            Lenght of the lattice.
            Only used if `lattice` is the name of a 1D Lattice.
        Lx, Ly : int
            Length of the lattice in x- and y-direction.
            Only used if `lattice` is the name of a 2D Lattice.
        bc_y : 'ladder' | 'cylinder'
            Boundary conditions in y-direction.
            Only used if `lattice` is the name of a 2D Lattice.
        """

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
        qvalue = get_parameter(model_params, 'phi', (1, 4), self.name)[1]
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        Lx = Lx_MUC * qvalue
        Ly = get_parameter(model_params, 'Ly', 4, self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = Square(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
        return lat

    def init_terms(self, model_params):

        t1 = get_parameter(model_params, 't1', 1., self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)
        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 4), self.name)
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        phi = 2 * np.pi * phi_p / phi_q

        u1, u2, dx = (0, 0, np.array([1, 0]))  # right
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext])
        self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.

        u1, u2, dx = (0, 0, np.array([0, 1]))  # up
        m = np.arange(0, phi_q * Lx_MUC)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j * phi * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.


class FermionicHofstadterModel(CouplingMPOModel):

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
        qvalue = get_parameter(model_params, 'phi', (1, 4), self.name)[1]
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        Lx = Lx_MUC * qvalue
        Ly = get_parameter(model_params, 'Ly', 4, self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = Square(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
        return lat

    def init_terms(self, model_params):
        t1 = get_parameter(model_params, 't1', 1., self.name, True)
        V = get_parameter(model_params, 'V', 10, self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)
        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 3), self.name)
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        phi = 2 * np.pi * phi_p / phi_q

        u1, u2, dx = (0, 0, np.array([1, 0]))  # right
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext])
        self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
        self.add_coupling(V, u1, 'N', u2, 'N', dx)

        u1, u2, dx = (0, 0, np.array([0, 1]))  # up
        m = np.arange(0, phi_q * Lx_MUC)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j * phi * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
        self.add_coupling(V, u1, 'N', u2, 'N', dx)


if __name__ == "__main__":

    t0 = time.time()

    model_params = dict(conserve='N', t1=1, filling=(int(1), int(9)), phi=(int(1), int(3)),
                        Lx_MUC=1, Ly=6, V=10,  # system params
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                        verbose=1, phi_ext=0)  # utility
    M = FermionicHofstadterModel(model_params)
    print("max MPO bond dimension = ", max(M.H_MPO.chi))

    print(time.time() - t0)
