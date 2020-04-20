# --- python imports
import matplotlib.pyplot as plt


#################################################################
# plot_MPS_unit_cell (plot the MPS unit cell for a given model) #
#################################################################


def plot_MPS_unit_cell(M):

    ax = plt.gca()
    # plotting
    M.lat.plot_basis(ax)
    M.lat.plot_bc_identified(ax)
    M.lat.plot_sites(ax)
    M.lat.plot_coupling(ax, linestyle='-', color='black')
    M.lat.plot_order(ax)
    # formatting
    plt.title('MPS unit cell')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')
    plt.show()
