# --- python imports
import numpy as np


#####################################################################
# scalar_observables (function for printing all scalar observables) #
#####################################################################


def scalar_observables(E, psi):
    print("Energy, E = ", E)
    print("Correlation length, xi = ", psi.correlation_length())
    print("von Neumann entanglement entropy, SvN =", psi.entanglement_entropy()[0])
    print("infinite Renyi entanglement entropy, Sinf = ", psi.entanglement_entropy(n=np.inf)[0])
    return


############################################################################################
# nonscalar_observables (function for recording the nonscalar observables listed in tools) #
############################################################################################


def nonscalar_observables(tools, data, psi, M, chiK_max, extra_dof, print_data=False):
    if 'ent_spec_real' in tools:
        ent_spec_real(data, psi, print_data)
    if 'ent_spec_mom' in tools:
        ent_spec_mom(data, psi, M, chiK_max, print_data)
    if 'density' in tools:
        expectation_value(data, 'density', psi, M, extra_dof, print_data)
    if 'corr_func' in tools:
        expectation_value(data, 'corr_func', psi, M, extra_dof, print_data)
    if 'corr_func_ext' in tools:
        expectation_value(data, 'corr_func_ext', psi, M, extra_dof, print_data)
    return


###############################################################################
# ent_spec_real (function for recording the real-space entanglement spectrum) #
###############################################################################


def ent_spec_real(data, psi, print_data):

    spectrum = psi.entanglement_spectrum(by_charge=True)

    print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=len(spectrum[0][0][0])))

    # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
    for bond in range(0, psi.L):
        for sector in range(0, len(spectrum[bond])):
            for i in range(0, len(spectrum[bond][sector][1])):
                data_line = "{charge:d}\t{bond:d}\t{spectrum:.15f}"\
                    .format(charge=spectrum[bond][sector][0][0],
                            bond=bond,
                            spectrum=spectrum[bond][sector][1][i])
                if print_data:
                    print(data_line)
                data['ent_spec_real'].write(data_line+"\n")

    return


##################################################################################
# ent_spec_mom (function for recording the momentum-space entanglement spectrum) #
##################################################################################


def ent_spec_mom(data, psi, M, chiK_max, print_data):

    (Un, W, q, ov, trunc_err) = \
        psi.compute_K(perm=M.lat, trunc_par={'chi_max': chiK_max}, canonicalize=1.e-6, verbose=0)

    if np.abs(np.abs(ov)-1) > 0.1:
        print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))
        print('Warning: State is not invariant under the permutation.')
    else:
        print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))

    print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=q.charges.shape[1]))

    # q.to_qflat()[i][0] --> q.to_qflat()[i][n] for different charge entries
    for i in range(len(W)):
        data_line = "{q:d}\t{K:.15f}\t{epsilon:.15f}"\
            .format(q=q.to_qflat()[i][0], K=np.angle(W[i])/np.pi, epsilon=-np.log(np.abs(W[i])))
        if print_data:
            print(data_line)
        data['ent_spec_mom'].write(data_line+"\n")

    return


###################################################################
# expectation_value (function for recording an expectation value) #
###################################################################

def expectation_value(data, name, psi, M, extra_dof, print_data):

    tot_numb_op = 'N' if not extra_dof else 'Ntot'

    if name == 'density':
        site_list = psi.expectation_value(tot_numb_op)
    elif name == 'corr_func':
        raw_site_list = psi.correlation_function(tot_numb_op, tot_numb_op, sites1=range(0, psi.L), sites2=[0])[:, 0]
        average_density = np.sum(psi.expectation_value(tot_numb_op))/psi.L
        site_list = [average_density - i for i in raw_site_list]
    elif name == 'corr_func_ext':
        raw_site_list = psi.correlation_function(tot_numb_op, tot_numb_op, sites1=range(0, 100*psi.L), sites2=[0])[:, 0]
        average_density = np.sum(psi.expectation_value(tot_numb_op)) / psi.L
        site_list = [average_density - i for i in raw_site_list]
        site_array = np.array(site_list)  # list type does not work in mps2lat_values_masked
        lattice_array = M.lat.mps2lat_values_masked(site_array)
    else:
        raise ValueError("Unknown name for the expectation_value function.")

    if not name.endswith('ext'):  # corr_func can be computed beyond psi.L, however then mps2lat function will not work
        lattice_array = M.lat.mps2lat_values(site_list)
        assert lattice_array.shape[:2] == M.lat.shape[:2]

    if print_data:
        print(lattice_array)
    if len(lattice_array.shape) == 2:  # one site per unit cell
        np.savetxt(data[name], lattice_array, delimiter='\t')
    elif len(lattice_array.shape) == 3:  # more than one site per unit cell
        for i in range(lattice_array.shape[-1]):
            np.savetxt(data[name], lattice_array[:, :, i], delimiter='\t')
    else:
        raise ValueError("Unexpected length of lattice_array in expectation_value function.")

    return
