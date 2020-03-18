# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_kappa_flow(threads, model, chi_max, use_pickle, make_pickle, **ham_params):

    fp.check_input_params("kappa_flow", threads, model, chi_max, ham_params, use_pickle, make_pickle)
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("kappa_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("kappa_flow", model, leaf)

    tools = ["corr_len_kappa_flow", "ent_spec_kappa_flow"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf)

    ####################################################################################################################

    t2, t2dash = ham_params['t2'], ham_params['t2dash']  # flow variables

    for kappa in np.linspace(ham_params['kappa_min'], ham_params['kappa_max'], ham_params['kappa_samp']):

        ham_params.update(t2=kappa*t2)
        ham_params.update(t2dash=kappa*t2dash)

        (E, psi, M) = fd.my_iDMRG_pickle("kappa_flow", model, chi_max, ham_params, use_pickle, make_pickle, run=True)

        #######################
        # corr_len_kappa_flow #
        #######################

        xi = psi.correlation_length()

        data_line = f"{kappa:.15f}\t{xi:.15f}"
        print(data_line)
        data['corr_len_kappa_flow'].write(data_line+"\n")

        #######################
        # ent_spec_kappa_flow #
        #######################

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        spectrum = psi.entanglement_spectrum(by_charge=True)

        bond = 0

        for sector in range(0, len(spectrum[bond])):
            for i in range(0, len(spectrum[bond][sector][1])):
                data_line = "{charge:d}\t{kappa:.15f}\t{spectrum:.15f}"\
                    .format(charge=spectrum[bond][sector][0][0],
                            kappa=kappa, spectrum=spectrum[bond][sector][1][i])
                print(data_line)
                data['ent_spec_kappa_flow'].write(data_line+"\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_kappa_flow(threads=1, model="FerHofHex1Hex5Orbital", chi_max=150, use_pickle=False, make_pickle=False,
                  t1=1, t5=-0.025, t5dash=0.1, kappa_min=0, kappa_max=1, kappa_samp=11, U=100, mu=0,
                  V=10, Vtype='Coulomb', Vrange=1, n=(1, 9), nphi=(1, 3), LxMUC=1, Ly=6, tag="")
