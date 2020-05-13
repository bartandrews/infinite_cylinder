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


def my_V_flow(path_flag, threads, model, chi_max, ham_params, use_pickle=False, make_pickle=False):

    path = "/home/bart/Desktop" if path_flag else ""  # specify the custom path
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("V_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("V_flow", path, model, chi_max, leaf)

    tools = ["corr_len_V_flow", "ent_spec_V_flow"]
    data = fp.prepare_output_files(tools, path, model, chi_max, leaf)

    ##################################################################################################################

    for V in np.linspace(ham_params['V_min'], ham_params['V_max'], ham_params['V_samp']):

        ham_params.update(V=V)
        (E, psi, M, _, _) = fd.my_iDMRG_pickle("V_flow", path, model, chi_max, ham_params, use_pickle, make_pickle,
                                               run=True)

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

    prog_args, stem_args, leaf_args = fa.parse_input_arguments("V_flow")

    my_V_flow(prog_args['path'], prog_args['threads'], stem_args['model'], stem_args['chi_max'], leaf_args,
              prog_args['use_pickle'], prog_args['make_pickle'])
