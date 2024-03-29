from tenpy.networks.mps import MPS
from tenpy.models.hubbard import FermiHubbardModel

from models.haldane import BosonicHaldaneModel, FermionicHaldaneModel
from models.hofstadter_old import BosonicHofstadterModel, FermionicHofstadterModel
from models.hofstadter_orbital import BosonicHofstadterOrbitalModel, FermionicHofstadterOrbitalModel

from models.hex_1_old import BosonicHex1Model, FermionicHex1Model
from models.hex_1_orbital import BosonicHex1OrbitalModel, FermionicHex1OrbitalModel
from models.tri_1 import BosonicTri1Model, FermionicTri1Model
from models.tri_2 import BosonicTri2Model, FermionicTri2Model
from models.hex_5 import BosonicHex5Model, FermionicHex5Model
from models.hex_1_hex_5_old import BosonicHex1Hex5Model, FermionicHex1Hex5Model
from models.hex_1_hex_5_orbital_old import BosonicHex1Hex5OrbitalModel, FermionicHex1Hex5OrbitalModel

from models.old.fermions_twist import FermionicTwistModel
from models.old.complete_twist import BosonicCompleteTwistModel, FermionicCompleteTwistModel

from models.examples.fermions_pi_flux import FermionicPiFluxModel
from models.examples.fermions_C3_haldane import FermionicC3HaldaneModel

from models.TBG.fermions_TBG1 import FermionicTBG1Model
from models.TBG.fermions_TBG2 import FermionicTBG2Model
from models.TBG.fermions_TBG3 import FermionicTBG3Model
from models.TBG.fermions_TBG4 import FermionicTBG4Model
from models.TBG.fermions_TBG5 import FermionicTBG5Model
from models.TBG.fermions_TBG6 import FermionicTBG6Model

from tenpy.algorithms import dmrg
# from tenpy.algorithms.mps_sweeps import OneSiteH, TwoSiteH

import sys
import random
import pickle


def file_name_stem(tool, model, lattice, initial_state, tile_unit, chi_max):

    if model not in ['Hubbard',
                     'BosonicHaldane', 'FermionicHaldane',
                     'BosonicHofstadter', 'FermionicHofstadter',
                     'BosonicHofstadterOrbital', 'FermionicHofstadterOrbital',
                     'BosonicHex1', 'FermionicHex1',
                     'BosonicHex1Orbital', 'FermionicHex1Orbital',
                     'BosonicTri1', 'FermionicTri1',
                     'BosonicTri2', 'FermionicTri2',
                     'BosonicHex5', 'FermionicHex5',
                     'BosonicHex1Hex5', 'FermionicHex1Hex5',
                     'BosonicHex1Hex5Orbital', 'FermionicHex1Hex5Orbital',
                     'FermionicTwist',
                     'BosonicCompleteTwist', 'FermionicCompleteTwist',
                     'FermionicPiFlux', 'FermionicC3Haldane',
                     'TBG1', 'TBG2', 'TBG3', 'TBG4', 'TBG5', 'TBG6']:
        sys.exit('Error: Unknown model.')

    stem = ("%s_%s_%s_%s_tile_%s_%s_chi_%s_"
            % (tool, model, lattice, initial_state, tile_unit[0], tile_unit[1], chi_max))

    return stem


def select_initial_psi(model, lattice, initial_state, tile_unit, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly):

    if lattice == "Square" or lattice == "Triangular":
        lat_basis = 1
    elif lattice == "Honeycomb":
        lat_basis = 2
    elif lattice == "TripartiteTriangular":
        lat_basis = 3
    elif lattice == "BipartiteSquare":
        lat_basis = 4
    elif lattice == "FiveBandLattice":
        lat_basis = 6
    elif lattice == "MagneticSquare":
        lat_basis = 4
    elif lattice == "MagneticHoneycomb":
        lat_basis = 8
    elif lattice == "MagneticTriangular":
        lat_basis = 4
    elif lattice == "MagneticTwist":
        lat_basis = 6
    else:
        sys.exit('Error: Unknown lattice.')

    Lx, Ly = model.lat.Ls[0], model.lat.Ls[1]

    if initial_state == 0:
        product_state = [tile_unit[0]] * model.lat.N_sites
    elif initial_state == 1:
        product_state = [tile_unit[1]] * model.lat.N_sites
    elif initial_state == 'random':
        product_state = []
        for i in range(model.lat.N_sites):
            product_state.append(random.choice([tile_unit[0], tile_unit[1]]))
    elif initial_state == 'neel':
        product_state = [tile_unit[0]]
        for i in range(0, lat_basis * Lx * Ly - 1):
            if i % 2 == 0:
                product_state.append(tile_unit[1])
            else:
                product_state.append(tile_unit[0])
    elif initial_state == 'third':
        product_state = [tile_unit[0]]
        for i in range(0, lat_basis * Lx * Ly - 1):
            if i % 3 == 0:
                product_state.append(tile_unit[1])
            else:
                product_state.append(tile_unit[0])
    elif initial_state == 'BosonicHofstadter' or initial_state == 'BosonicTri1':
        if (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 4):  # nphi = 1/4 (Ly=4)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 6):  # nphi = 1/4 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 4):  # nphi = 1/5 (Ly=4)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 6):  # nphi = 1/5 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx, Ly) == (1, 4):  # nphi = 1/6 (Ly=4)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx, Ly) == (1, 6):  # nphi = 1/6 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            sys.exit('Error: Unknown initial_state configuration.')
    elif initial_state == 'BosonicHofstadterOrbital':
        product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0,
                         '1_x 0_y', 0, 0, 0, 0, 0, 0, 0]
    elif initial_state == 'FermionicHofstadter' or initial_state == 'FermionicTri1':
        if (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx, Ly) == (1, 6):  # nphi = 1/3 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx, Ly) == (1, 9):  # nphi = 1/3 (Ly=9)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 6):  # nphi = 1/4 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 9):  # nphi = 1/4 (Ly=9)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 6):  # nphi = 1/5 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 9):  # nphi = 1/5 (Ly=9)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            sys.exit('Error: Unknown initial_state configuration.')
    elif initial_state == 'FermionicHofstadterOrbital':
        product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0,
                         'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0]
    elif initial_state == 'BosonicHex1' or initial_state == 'BosonicHex5' or initial_state == 'BosonicHex1Hex5':
        if (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 4):  # nphi = 1/4 (Ly=4)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 4):  # nphi = 1/4 (Ly=4) CI
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 6):  # nphi = 1/4 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 4):  # nphi = 1/5 (Ly=4)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 5) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 4):  # nphi = 1/5 (Ly=4) CI
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 6):  # nphi = 1/5 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx, Ly) == (1, 4):  # nphi = 1/6 (Ly=4)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 6) and (pvalue, qvalue) == (1, 6) and (Lx, Ly) == (1, 4):  # nphi = 1/6 (Ly=4) CI
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx, Ly) == (1, 6):  # nphi = 1/6 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            sys.exit('Error: Unknown initial_state configuration.')
    elif initial_state == 'BosonicHex1Orbital':
        product_state = ['1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         '1_x 0_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    elif initial_state == 'FermionicHex1' or initial_state == 'FermionicHex5' or initial_state == 'FermionicHex1Hex5':
        if (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx, Ly) == (1, 6):  # nphi = 1/3 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 3) and (pvalue, qvalue) == (1, 3) and (Lx, Ly) == (1, 6):  # nphi = 1/3 (Ly=6) CI
            product_state = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx, Ly) == (1, 9):  # nphi = 1/3 (Ly=9)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 6):  # nphi = 1/4 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 6):  # nphi = 1/4 (Ly=6) CI
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 9):  # nphi = 1/4 (Ly=9)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 6):  # nphi = 1/5 (Ly=6)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 5) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 6):  # nphi = 1/5 (Ly=6) CI
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 9):  # nphi = 1/5 (Ly=9)
            product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            sys.exit('Error: Unknown initial_state configuration.')
    elif initial_state == 'FermionicHex1Orbital':
        product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    elif initial_state == 'BosonicHex1Hex5Orbital':
        if (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 4):  # nphi = 1/4 (Ly=4)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 4):  # nphi = 1/4 (Ly=4) CI
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 8) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 6):  # nphi = 1/4 (Ly=6)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 4):  # nphi = 1/5 (Ly=4)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 5) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 4):  # nphi = 1/5 (Ly=4) CI
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 10) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 6):  # nphi = 1/5 (Ly=6)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx, Ly) == (1, 4):  # nphi = 1/6 (Ly=4)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 6) and (pvalue, qvalue) == (1, 6) and (Lx, Ly) == (1, 4):  # nphi = 1/6 (Ly=4) CI
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 6) and (Lx, Ly) == (1, 6):  # nphi = 1/6 (Ly=6)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            sys.exit('Error: Unknown initial_state configuration.')
    elif initial_state == 'FermionicHex1Hex5Orbital':
        if (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx, Ly) == (1, 6):  # nphi = 1/3 (Ly=6)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 3) and (pvalue, qvalue) == (1, 3) and (Lx, Ly) == (1, 6):  # nphi = 1/3 (Ly=6) CI
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 9) and (pvalue, qvalue) == (1, 3) and (Lx, Ly) == (1, 9):  # nphi = 1/3 (Ly=9)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 6):  # nphi = 1/4 (Ly=6)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 4) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 6):  # nphi = 1/4 (Ly=6) CI
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 12) and (pvalue, qvalue) == (1, 4) and (Lx, Ly) == (1, 9):  # nphi = 1/4 (Ly=9)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 6):  # nphi = 1/5 (Ly=6)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 5) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 6):  # nphi = 1/5 (Ly=6) CI
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0]
        elif (nnvalue, ndvalue) == (1, 15) and (pvalue, qvalue) == (1, 5) and (Lx, Ly) == (1, 9):  # nphi = 1/5 (Ly=9)
            product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        else:
            sys.exit('Error: Unknown initial_state configuration.')
    else:
        sys.exit('Error: Unknown initial_state.')

    return product_state


def define_iDMRG_model(model, lattice, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, phi_ext=0, invt2dash=0):

    if model == 'Hubbard':
        model_params = dict(cons_N='N', cons_Sz='Sz', t=t1, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermiHubbardModel(model_params)

    elif model == 'BosonicHaldane':
        model_params = dict(conserve='N', t=t1, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0, phi_ext=phi_ext)
        M = BosonicHaldaneModel(model_params)

    elif model == 'FermionicHaldane':
        model_params = dict(conserve='N', t=t1, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0, phi_ext=phi_ext)
        M = FermionicHaldaneModel(model_params)

    elif model == 'BosonicHofstadter':
        model_params = dict(conserve='N', t=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHofstadterModel(model_params)

    elif model == 'FermionicHofstadter':
        model_params = dict(conserve='N', t=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, V=10,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHofstadterModel(model_params)

    elif model == 'BosonicHofstadterOrbital':
        model_params = dict(conserve='N', t=t1, filling=(1, 8), phi=(1, 4), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHofstadterOrbitalModel(model_params)

    elif model == 'FermionicHofstadterOrbital':
        model_params = dict(conserve='N', t=t1, filling=(1, 9), phi=(1, 3), Lx=Lx, Ly=Ly, V=10,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHofstadterOrbitalModel(model_params)

    ####################################################################################################################

    elif model == 'BosonicHex1':
        model_params = dict(conserve='N', t=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHex1Model(model_params)

    elif model == 'FermionicHex1':
        model_params = dict(conserve='N', t=t1, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHex1Model(model_params)

    elif model == 'BosonicHex1Orbital':
        model_params = dict(conserve='N', t=t1, filling=(1, 8), phi=(1, 4), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHex1OrbitalModel(model_params)

    elif model == 'FermionicHex1Orbital':
        model_params = dict(conserve='N', t=t1, filling=(1, 9), phi=(1, 3), Lx=Lx, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHex1OrbitalModel(model_params)

    elif model == 'BosonicTri1':
        model_params = dict(conserve='N', t=t1, filling=(1, 8), phi=(1, 4), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicTri1Model(model_params)

    elif model == 'FermionicTri1':
        model_params = dict(conserve='N', t=t1, filling=(1, 9), phi=(1, 3), Lx=Lx, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicTri1Model(model_params)

    elif model == 'BosonicTri2':
        model_params = dict(conserve='N', t=t1, filling=(1, 8), phi=(1, 4), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicTri2Model(model_params)

    elif model == 'FermionicTri2':
        model_params = dict(conserve='N', t=t1, filling=(1, 9), phi=(1, 3), Lx=Lx, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicTri2Model(model_params)

    elif model == 'BosonicHex5':
        model_params = dict(conserve='N', t=t2, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHex5Model(model_params)

    elif model == 'FermionicHex5':
        model_params = dict(conserve='N', t=t2, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHex5Model(model_params)

    elif model == 'BosonicHex1Hex5':
        model_params = dict(conserve='N', t1=t1, t2=t2, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHex1Hex5Model(model_params)

    elif model == 'FermionicHex1Hex5':
        model_params = dict(conserve='N', t1=t1, t2=t2, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = FermionicHex1Hex5Model(model_params)

    elif model == 'BosonicHex1Hex5Orbital':
        model_params = dict(conserve='N', t1=t1, t2=t2, t2dash=t2dash, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, Nmax=1,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext)  # utility
        M = BosonicHex1Hex5OrbitalModel(model_params)

    elif model == 'FermionicHex1Hex5Orbital':
        model_params = dict(conserve='N', t1=t1, t2=t2, t2dash=t2dash, filling=(int(nnvalue), int(ndvalue)), phi=(int(pvalue), int(qvalue)), Lx=Lx, Ly=Ly, V=V,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=phi_ext, invt2dash=invt2dash)  # utility
        M = FermionicHex1Hex5OrbitalModel(model_params)

    ####################################################################################################################

    elif model == 'FermionicTwist':
        model_params = dict(conserve='N', t=t1, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicTwistModel(model_params)

    elif model == 'BosonicCompleteTwist':
        model_params = dict(conserve='N', t=t1, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = BosonicCompleteTwistModel(model_params)

    elif model == 'FermionicCompleteTwist':
        model_params = dict(conserve='N', t=t1, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicCompleteTwistModel(model_params)

    elif model == 'FermionicPiFlux':
        model_params = dict(conserve='N', t=t1, mu=mu, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicPiFluxModel(model_params)

    elif model == 'FermionicC3Haldane':
        model_params = dict(conserve='N', t=t1, mu=mu, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1, phi_ext=phi_ext)
        M = FermionicC3HaldaneModel(model_params)

    elif model == 'TBG1':
        model_params = dict(cons_N='N', cons_Sz='Sz', t=t1, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG1Model(model_params)

    elif model == 'TBG2':
        model_params = dict(conserve='N', t=t1, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG2Model(model_params)

    elif model == 'TBG3':
        model_params = dict(cons_N='N', cons_Sz='Sz', t=t1, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG3Model(model_params)

    elif model == 'TBG4':
        model_params = dict(conserve='N', t=t1, U=U, mu=mu, V=V, lattice=lattice, Lx=Lx, Ly=Ly, verbose=1)
        M = FermionicTBG4Model(model_params)

    elif model == 'TBG6':
        model_params = dict(conserve='N', t=t1, U=U, mu=mu, V=V, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=0)
        M = FermionicTBG6Model(model_params)

    return M


def define_iDMRG_spin_model(model, lattice, J, Js, Jv, Lx, Ly, phi_ext=0):

    if model == 'TBG5':
        model_params = dict(conserve='Sz', J=J, Js=Js, Jv=Jv, lattice=lattice, bc_MPS='infinite',
                            order='default', Lx=Lx, Ly=Ly, bc_y='cylinder', verbose=1, phi_ext=phi_ext)
        M = FermionicTBG5Model(model_params)

    return M


def define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, phi_ext=0, invt2dash=0):

    M = define_iDMRG_model(model, lattice, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, phi_ext, invt2dash)

    product_state = select_initial_psi(M, lattice, initial_state, tile_unit, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    # refined

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True,
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
        'max_E_err': 1.e-10,
        'max_S_err': 1.e-6,
        # 'norm_tol': 1.e-6,
        # 'norm_tol_iter': 1000,
        # 'max_sweeps': 150,
        'verbose': 1,
        'N_sweeps_check': 10
    }

    # unrefined

    # dmrg_params = {
    #     'mixer': True,  # setting this to True helps to escape local minima
    #     'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
    #     'trunc_params': {
    #         'chi_max': chi_max,
    #         'svd_min': 1.e-10
    #     },
    #     # 'lanczos_params': {
    #     #     'reortho': True,
    #     #     'N_cache': 40
    #     # },
    #     # 'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
    #     'max_E_err': 1.e-8,
    #     'max_S_err': 1.e-6,
    #     # 'norm_tol': 1.e-6,
    #     # 'norm_tol_iter': 1000,
    #     'max_sweeps': 150,
    #     'verbose': 1,
    #     # 'N_sweeps_check': 10
    # }

    # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
    engine = dmrg.EngineCombine(psi, M, dmrg_params)

    return engine


def define_iDMRG_spin_engine(model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly, phi_ext=0):

    M = define_iDMRG_spin_model(model, lattice, J, Js, Jv, Lx, Ly, phi_ext)

    product_state = select_initial_psi(M, lattice, initial_state, tile_unit)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True,
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
        'max_E_err': 1.e-10,
        'max_S_err': 1.e-6,
        # 'norm_tol': 1.e-6,
        # 'norm_tol_iter': 1000,
        # 'max_sweeps': 150,
        'verbose': 1,
        'N_sweeps_check': 10
    }

    engine = dmrg.OneSiteDMRGEngine(psi, M, dmrg_params)

    return engine


def define_iDMRG_engine_pickle(flow, model, lattice, initial_state, tile_unit, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly,
                               use_pickle=False, make_pickle=False, phi_ext=0, invt2dash=0):

    global qval
    qval = qvalue

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("engine", model, lattice, initial_state, tile_unit, chi_max)
        pickle_leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_n_%s_%s_nphi_%s_%s_Lx_%s_Ly_%s_phi_%s.pkl" % (t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:

        with open(pickle_file, 'rb') as file1:
            engine = pickle.load(file1)

    else:

        engine = define_iDMRG_engine(model, lattice, initial_state, tile_unit, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, phi_ext,
                                     invt2dash)

        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                pickle.dump(engine, file2)

    return engine


def define_iDMRG_spin_engine_pickle(flow, model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly,
                                    use_pickle=False, make_pickle=False, phi_ext=0):

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("engine", model, lattice, initial_state, tile_unit, chi_max)
        pickle_leaf = ("J_%s_Js_%s_Jv_%s_Lx_%s_Ly_%s_phi_%s.pkl" % (J, Js, Jv, Lx, Ly, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:

        with open(pickle_file, 'rb') as file1:
            engine = pickle.load(file1)

    else:

        engine = define_iDMRG_spin_engine(model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly, phi_ext)

        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                pickle.dump(engine, file2)

    return engine


def run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t1, t2, t2dash, U, mu, V, Lx, Ly, nnvalue, ndvalue, pvalue, qvalue, phi_ext=0, invt2dash=0):

    M = define_iDMRG_model(model, lattice, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, phi_ext, invt2dash)

    product_state = select_initial_psi(M, lattice, initial_state, tile_unit, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly)

    print(product_state)  # NB: two sites per basis for honeycomb crystal

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    # refined

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
        'max_E_err': 1.e-10,
        'max_S_err': 1.e-6,
        # 'max_sweeps': 150,
        'verbose': 1,
        'N_sweeps_check': 10
    }

    # unrefined

    # dmrg_params = {
    #     'mixer': True,  # setting this to True helps to escape local minima
    #     'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
    #     'trunc_params': {
    #         'chi_max': chi_max,
    #         'svd_min': 1.e-10
    #     },
    #     # 'lanczos_params': {
    #     #     'reortho': True,
    #     #     'N_cache': 40
    #     # },
    #     # 'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
    #     'max_E_err': 1.e-8,
    #     'max_S_err': 1.e-6,
    #     # 'norm_tol': 1.e-6,
    #     # 'norm_tol_iter': 1000,
    #     'max_sweeps': 150,
    #     'verbose': 1,
    #     # 'N_sweeps_check': 10
    # }

    info = dmrg.run(psi, M, dmrg_params)

    E = info['E']

    return E, psi, M


# def run_iDMRG_spin(model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly, phi_ext=0):
#
#     M = define_iDMRG_spin_model(model, lattice, J, Js, Jv, Lx, Ly, phi_ext)
#
#     product_state = select_initial_psi(M, lattice, initial_state, tile_unit)
#
#     print(product_state)  # NB: two sites per basis for honeycomb crystal
#
#     psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)
#
#     dmrg_params = {
#         'mixer': True,  # setting this to True helps to escape local minima
#         'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
#         'trunc_params': {
#             # 'chi_max': chi_max,
#             'svd_min': 1.e-10
#         },
#         # 'lanczos_params': {
#         #     'reortho': True
#         #     'N_cache': 40
#         # },
#         'chi_list': {0: 9, 10: 49, 20: 100, 40: chi_max},
#         'max_E_err': 1.e-10,
#         'max_S_err': 1.e-6,
#         # 'max_sweeps': 150,
#         'verbose': 1,
#         'N_sweeps_check': 10
#     }
#
#     info = dmrg.run(psi, M, dmrg_params)
#
#     E = info['E']
#
#     return E, psi, M


def run_iDMRG_pickle(flow, model, lattice, initial_state, tile_unit, chi_max, t1, t2, t2dash, U, mu, V, Lx, Ly, nnvalue, ndvalue, pvalue, qvalue,
                     use_pickle=False, make_pickle=False, phi_ext=0, invt2dash=0):

    # update q for corresponding Lattice class
    global qval
    qval = qvalue

    if use_pickle or make_pickle:
        pickle_stem = file_name_stem("E_psi_M", model, lattice, initial_state, tile_unit, chi_max)
        pickle_leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_n_%s_%s_nphi_%s_%s_phi_%s.pkl" % (t1, t2, t2dash, U, mu, V, Lx, Ly, nnvalue, ndvalue, pvalue, qvalue, phi_ext))
        pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf

    if use_pickle:

        with open(pickle_file, 'rb') as file1:
            [E, psi, M] = pickle.load(file1)

    else:

        (E, psi, M) = run_iDMRG(model, lattice, initial_state, tile_unit, chi_max,
                                t1, t2, t2dash, U, mu, V, Lx, Ly, nnvalue, ndvalue, pvalue, qvalue, phi_ext, invt2dash)

        if make_pickle:
            with open(pickle_file, 'wb') as file2:
                pickle.dump([E, psi, M], file2)

    return E, psi, M


# def run_iDMRG_spin_pickle(flow, model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly,
#                           use_pickle=False, make_pickle=False, phi_ext=0):
#
#     if use_pickle or make_pickle:
#         pickle_stem = file_name_stem("E_psi_M", model, lattice, initial_state, tile_unit, chi_max)
#         pickle_leaf = ("J_%s_Js_%s_Jv_%s_Lx_%s_Ly_%s_phi_%s.pkl" % (J, Js, Jv, Lx, Ly, phi_ext))
#         pickle_file = "pickles/" + flow + "/" + pickle_stem + pickle_leaf
#
#     if use_pickle:
#
#         with open(pickle_file, 'rb') as file1:
#             [E, psi, M] = pickle.load(file1)
#
#     else:
#
#         (E, psi, M) = run_iDMRG_spin(model, lattice, initial_state, tile_unit, chi_max, J, Js, Jv, Lx, Ly, phi_ext)
#
#         if make_pickle:
#             with open(pickle_file, 'wb') as file2:
#                 pickle.dump([E, psi, M], file2)
#
#     return E, psi, M

