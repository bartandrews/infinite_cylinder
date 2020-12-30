# --- python imports
import matplotlib.pyplot as plt


#################################################################
# plot_MPS_unit_cell (plot the MPS unit cell for a given model) #
#################################################################


def plot_MPS_unit_cell(M):

    ax = plt.gca()
    # plotting
    # M.lat.plot_basis(ax)
    # M.lat.plot_bc_identified(ax)
    M.lat.plot_sites(ax)
    # M.lat.plot_coupling(ax, linestyle='-', color='black')
    # M.coupling_terms['t1'].plot_coupling_terms(ax, M.lat, text='t1={strength:.2g}')
    # M.coupling_terms['t2'].plot_coupling_terms(ax, M.lat, text='t2={strength:.2g}')
    # M.coupling_terms['t3'].plot_coupling_terms(ax, M.lat, text='t3={strength:.2g}')
    M.coupling_terms['t4'].plot_coupling_terms(ax, M.lat, text='t4={strength:.2g}')
    # M.lat.plot_order(ax)
    # formatting
    plt.title('MPS unit cell')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')
    plt.show()
