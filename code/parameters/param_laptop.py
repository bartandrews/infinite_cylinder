# configuration parameters
model = 'TBG4'
lattice = 'FiveBandLattice'
initial_state = 'custom'

if model == 'BosonicHaldane':
    tile_unit = ['0', '1']
elif model == 'BosonicHaldane2':
    tile_unit = ['0', '1']
elif model == 'FermionicHaldane':
    tile_unit = ['empty', 'full']
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

# chi_max for DMRG
chi_max = 500
# chi max for compute_K
chi_max_K = 1000

# Hamiltonian parameters (U=0 for FermionicHaldane)
t, mu, V = -1, 0, 1

if model in ['BosonicHaldane', 'BosonicHaldane2', 'FermionicHaldane']:
    U = 0
elif model in ['Hubbard', 'TBG1', 'TBG2']:
    U = 1
elif model in ['TBG3', 'TBG4']:
    mu, V, U = 0, 0, 0

# unit cell
Lx, Ly = 1, 3

# pickle capability
use_pickle = False
make_pickle = False
