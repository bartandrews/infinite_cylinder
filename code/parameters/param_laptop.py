# configuration parameters
model = 'FermionicHaldane'
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

chi_max = 10
# Hamiltonian parameters (U=0 for FermionicHaldane)
t, mu, V = -1, 0, 1

if model in ['BosonicHaldane', 'FermionicHaldane']:
    U = 0
elif model in ['Hubbard', 'TBG1', 'TBG2']:
    U = 0

# unit cell
Lx, Ly = 1, 3
