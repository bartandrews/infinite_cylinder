import numpy as np
import sys
import time

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite, GroupedSite
from tenpy.models.lattice import Honeycomb


def get_product_state(model, nnvalue, ndvalue, qvalue, Lx_MUC, Ly, filling_scale_factor=1, orbital_preference=None):

    if "Orbital" in model:
        numb_particles = 2 * int(qvalue) * int(Lx_MUC) * int(Ly) * int(filling_scale_factor) * int(nnvalue) / int(ndvalue)
    else:
        numb_particles = int(qvalue) * int(Lx_MUC) * int(Ly) * int(filling_scale_factor) * int(nnvalue) / int(ndvalue)

    if not numb_particles.is_integer():
        sys.exit("Error: Cannot fit an integer number of particles into lattice geometry.")
    else:
        numb_particles = int(numb_particles)

    if "Hofstadter" in model:
        system_size = int(qvalue * Lx_MUC * Ly)
    elif "Hex" in model:
        system_size = int(2 * qvalue * Lx_MUC * Ly)
    numb_sites_per_particle = int(system_size / numb_particles)

    if "Orbital" in model:
        if orbital_preference in ['polarized_x', None]:
            if "Bosonic" in model:
                lattice_site = ['1_x 0_y']
            elif "Fermionic" in model:
                lattice_site = ['full_x empty_y']
            state = (lattice_site + [0] * (numb_sites_per_particle - 1)) * int((system_size / numb_sites_per_particle))
        elif orbital_preference == 'polarized_y':
            if "Bosonic" in model:
                lattice_site = ['0_x 1_y']
            elif "Fermionic" in model:
                lattice_site = ['empty_x full_y']
            state = (lattice_site + [0] * (numb_sites_per_particle - 1)) * int((system_size / numb_sites_per_particle))
        elif orbital_preference == 'unpolarized':
            if "Bosonic" in model:
                lattice_site_1 = ['1_x 0_y']
                lattice_site_2 = ['0_x 1_y']
            elif "Fermionic" in model:
                lattice_site_1 = ['full_x empty_y']
                lattice_site_2 = ['empty_x full_y']
            state = (lattice_site_1 + [0] * (numb_sites_per_particle - 1)
                     + lattice_site_2 + [0] * (numb_sites_per_particle - 1)) \
                    * int((system_size / (2 * numb_sites_per_particle)))
        elif orbital_preference == 'filled':
            if "Bosonic" in model:
                lattice_site = ['1_x 1_y']
            elif "Fermionic" in model:
                lattice_site = ['full_x full_y']
            state = (lattice_site + [0] * (2 * numb_sites_per_particle - 1)) \
                    * int((system_size / (2 * numb_sites_per_particle)))
        else:
            sys.exit("Error: Unknown orbital_preference parameter.")
    else:
        state = ([1] + [0] * (numb_sites_per_particle - 1)) * int((system_size / numb_sites_per_particle))

    print("initial state = ", state)
    print("number of particles = ", numb_particles)
    print("number of lattice sites = ", len(state))
    if "Orbital" in model:
        print("number of orbital sites = ", 2*len(state))

    return state


if __name__ == "__main__":

    state = get_product_state(model="FermionicHex1Hex5Orbital", nnvalue=1, ndvalue=15, qvalue=3, Lx_MUC=1, Ly=5,
                              filling_scale_factor=1, orbital_preference='polarized_x')

