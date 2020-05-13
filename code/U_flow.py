# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd
import functions.func_args as fa


def my_U_flow(path_flag, threads, model, chi_max, ham_params):

    path = "/home/bart/Desktop" if path_flag else ""  # specify the custom path
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("U_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("U_flow", path, model, chi_max, leaf)

    tools = ["corr_len_U_flow", "double_occ_U_flow"]
    data = fp.prepare_output_files(tools, path, model, chi_max, leaf)

    ####################################################################################################################

    for U in np.linspace(ham_params['U_min'], ham_params['U_max'], ham_params['U_samp']):

        ham_params.update(U=U)
        (E, psi, M, _, _) = fd.my_iDMRG_pickle("U_flow", path, model, chi_max, ham_params, run=True)

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

    prog_args, stem_args, leaf_args = fa.parse_input_arguments("U_flow")

    my_U_flow(prog_args['path'], prog_args['threads'], stem_args['model'], stem_args['chi_max'], leaf_args)
