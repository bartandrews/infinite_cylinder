# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_V_flow(threads, model, chi_max, use_pickle=False, make_pickle=False, **ham_params):

    fp.check_input_params("V_flow", threads, model, chi_max, ham_params, use_pickle, make_pickle)
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("V_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("V_flow", model, leaf)

    tools = ["corr_len_V_flow", "ent_spec_V_flow"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf)

    ##################################################################################################################

    for V in np.linspace(ham_params['V_min'], ham_params['V_max'], ham_params['V_samp']):

        ham_params.update(V=V)
        (E, psi, M) = fd.my_iDMRG_pickle("V_flow", model, chi_max, ham_params, use_pickle, make_pickle, run=True)

        ###################
        # corr_len_V_flow #
        ###################

        xi = psi.correlation_length()

        data_line = f"{V:.15f}\t{xi:.15f}"
        print(data_line)
        data['corr_len_V_flow'].write(data_line+"\n")

        ###################
        # ent_spec_V_flow #
        ###################

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        spectrum = psi.entanglement_spectrum(by_charge=True)

        bond = 0

        for sector in range(0, len(spectrum[bond])):
            for i in range(0, len(spectrum[bond][sector][1])):
                data_line = "{charge:d}\t{V:.15f}\t{spectrum:.15f}"\
                    .format(charge=spectrum[bond][sector][0][0],
                            V=V, spectrum=spectrum[bond][sector][1][i])
                print(data_line)
                data['ent_spec_V_flow'].write(data_line+"\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_V_flow(threads=1, model="FerHofHex1Hex5Orbital", chi_max=50, use_pickle=False, make_pickle=False,
              t1=1, t5=-0.025, t5dash=0.1, U=0, mu=0, V_min=0, V_max=10, V_samp=100, Vtype='Coulomb', Vrange=1,
              nn=1, nd=9, p=1, q=3, LxMUC=1, Ly=6, tag="")
