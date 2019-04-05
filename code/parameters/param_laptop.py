# configuration parameters
model = 'TBG3'
lattice = 'Honeycomb'
initial_state = 'neel'

if model == 'BosonicHaldane':
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

# chi_max for DMRG
chi_max = 500
# chi max for compute_K
chi_max_K = 100

# Hamiltonian parameters (U=0 for FermionicHaldane)
t, mu, V = -1, 0, 1

if model in ['BosonicHaldane', 'FermionicHaldane']:
    U = 0
elif model in ['Hubbard', 'TBG1', 'TBG2', 'TBG3']:
    U = 1

# unit cell
Lx, Ly = 1, 4

# pickle capability
use_pickle = False
make_pickle = False
