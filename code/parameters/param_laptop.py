# name tag
tag = 'boson_test'

# configuration parameters
model = 'BosonicHofstadter'
lattice = 'MagneticSquare'
initial_state = 'BosonicHofstadter'

if model == 'Hubbard':
    tile_unit = ['down', 'up']
elif model == 'BosonicHaldane' or model == 'FermionicHaldane':
    tile_unit = [0, 1]
elif model == 'BosonicHofstadter' or model == 'FermionicHofstadter':
    tile_unit = [0, 1]
elif model == 'BosonicHofstadterOrbital' or model == 'FermionicHofstadterOrbital':
    tile_unit = ['full_x empty_y', 'full_x empty_y']
########################################################################################################################
elif model == 'BosonicHex1' or model == 'FermionicHex1':
    tile_unit = [0, 1]
elif model == 'BosonicHex1Orbital' or model == 'FermionicHex1Orbital':
    tile_unit = ['full_x empty_y', 'full_x empty_y']
elif model == 'BosonicTri1' or model == 'FermionicTri1':
    tile_unit = [0, 1]
elif model == 'BosonicTri2' or model == 'FermionicTri2':
    tile_unit = [0, 1]
elif model == 'BosonicHex5' or model == 'FermionicHex5':
    tile_unit = [0, 1]
elif model == 'BosonicHex1Hex5' or model == 'FermionicHex1Hex5':
    tile_unit = [0, 1]
elif model == 'BosonicHex1Hex5Orbital' or model == 'FermionicHex1Hex5Orbital':
    tile_unit = ['empty_x empty_y', 'full_x full_y']
########################################################################################################################
elif model == 'FermionicTwist':
    tile_unit = [0, 1]
elif model == 'BosonicCompleteTwist' or model == 'FermionicCompleteTwist':
    tile_unit = ['empty_x empty_y', 'full_x full_y']
elif model == 'FermionicPiFlux':
    tile_unit = ['empty', 'full']
elif model == 'FermionicC3Haldane':
    tile_unit = ['full_A empty_B', 'empty_A full_B', 'full_A empty_B']
elif model == 'TBG1':
    tile_unit = ['down_px up_py', 'down_px up_py']
elif model == 'TBG2':
    tile_unit = ['empty_px full_py', 'empty_px full_py']
elif model == 'TBG3':
    tile_unit = ['down_px empty_py', 'empty_px up_py']
elif model == 'TBG4':
    tile_unit = ['full_pz empty_pp empty_pm', 'full', 'full', 'full_pz empty_pp empty_pm', 'full', 'full']
elif model == 'TBG5':
    tile_unit = ['up_spin down_valley', 'up_spin up_valley']
elif model == 'TBG6':
    tile_unit = ['full_px empty_py full_z', 'empty_px full_py empty_z']

# chi_max for DMRG
chi_max = 50
# chi max for compute_K
chi_max_K = 50

# Hamiltonian parameters (U=0 for FermionicHaldane)
t, mu, U, V = -1, 0, 0, 0

if model in ['FermionicHaldane']:
    U = 0
elif model in ['Hubbard', 'TBG1', 'TBG2']:
    U = 1
elif model in ['TBG3', 'TBG4']:
    mu, V, U = 0, 0, 0
elif model in ['TBG5']:
    J, Js, Jv = 1, 0.1, 8

# unit cell
Lx, Ly = 1, 4

# pickle capability
use_pickle = False
make_pickle = False
