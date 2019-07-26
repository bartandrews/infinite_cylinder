# configuration parameters
model = 'FermionicTwist'
lattice = 'MagneticTwist'
initial_state = 'third'

if model == 'BosonicHaldane':
    tile_unit = ['0', '1']
elif model == 'BosonicHaldane2':
    tile_unit = ['0', '1']
elif model == 'FermionicHaldane':
    tile_unit = ['empty', 'full']
elif model == 'FermionicHofstadter':
    tile_unit = ['empty', 'full']
elif model == 'FermionicHex1':
    tile_unit = ['empty', 'full']
elif model == 'FermionicTwist':
    tile_unit = ['empty', 'full']
elif model == 'FermionicPiFlux':
    tile_unit = ['empty', 'full']
elif model == 'FermionicC3Haldane':
    tile_unit = ['full_A empty_B', 'empty_A full_B', 'full_A empty_B']
elif model == 'Hubbard':
    tile_unit = ['down', 'up']
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
chi_max = 500
# chi max for compute_K
chi_max_K = 500

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
Lx, Ly = 3, 1

# pickle capability
use_pickle = False
make_pickle = False
