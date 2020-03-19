# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_U_flow(threads, model, chi_max, use_pickle=False, make_pickle=False, **ham_params):

    fp.check_input_params("U_flow", threads, model, chi_max, ham_params, use_pickle, make_pickle)
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("U_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("U_flow", model, leaf)

    tools = ["corr_len_U_flow", "double_occ_U_flow"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf)

    ####################################################################################################################

    for U in np.linspace(ham_params['U_min'], ham_params['U_max'], ham_params['U_samp']):

        ham_params.update(U=U)
        (E, psi, M) = fd.my_iDMRG_pickle("U_flow", model, chi_max, ham_params, use_pickle, make_pickle, run=True)

        ###################
        # corr_len_U_flow #
        ###################

        xi = psi.correlation_length()

        data_line = f"{U:.15f}\t{xi:.15f}"
        print(data_line)
        data['corr_len_U_flow'].write(data_line+"\n")

        #####################
        # double_occ_U_flow #
        #####################

        term1 = psi.expectation_value('Ntot Ntot')
        term2 = -2*psi.expectation_value('Ntot')

        term1_2 = term1 + term2

        nd_summand = [x + 1 for x in term1_2]
        nd = np.average(nd_summand)

        data_line = f"{U:.15f}\t{nd:.15f}"
        print(data_line)
        data['double_occ_U_flow'].write(data_line+"\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_U_flow(threads=1, model="FerHofHex1Hex5Orbital", chi_max=50, use_pickle=False, make_pickle=False,
              t1=1, t5=-0.025, t5dash=0.1, U_min=0, U_max=10, U_samp=100, mu=0, V=0, Vtype='Coulomb', Vrange=1,
              nn=1, nd=9, p=1, q=3, LxMUC=1, Ly=6, tag="")
