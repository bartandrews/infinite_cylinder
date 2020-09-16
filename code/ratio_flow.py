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


def my_ratio_flow(path_flag, threads, model, chi_max, ham_params):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("ratio_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("ratio_flow", path, model, chi_max, leaf)

    tools = ["ent_spec_ratio_flow", "ent_ratio_flow", "corr_len_ratio_flow", "ent_corr_len_ratio"]
    data = fp.prepare_output_files(tools, path, model, chi_max, leaf)

    ##################################################################################################################

    for r in np.linspace(ham_params['r_min'], ham_params['r_max'], ham_params['r_samp']):

        ham_params.update(r=r)
        state_data = fd.my_iDMRG_pickle("ratio_flow", path, model, chi_max, ham_params, run=True)
        psi = state_data['psi']

        #######################
        # ent_spec_ratio_flow #
        #######################

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        spectrum = psi.entanglement_spectrum(by_charge=True)

        bond = 0

        for sector in range(0, len(spectrum[bond])):
            for i in range(0, len(spectrum[bond][sector][1])):
                data_line = "{charge:d}\t{ratio:.15f}\t{spectrum:.15f}"\
                    .format(charge=spectrum[bond][sector][0][0],
                            ratio=r, spectrum=spectrum[bond][sector][1][i])
                print(data_line)
                data['ent_spec_ratio_flow'].write(data_line+"\n")

        ##################
        # ent_ratio_flow #
        ##################

        SvN = psi.entanglement_entropy()[0]

        data_line = f"{r:.15f}\t{SvN:.15f}"
        print(data_line)
        data['ent_ratio_flow'].write(data_line + "\n")

        #######################
        # corr_len_ratio_flow #
        #######################

        xi = psi.correlation_length()

        data_line = f"{r:.15f}\t{xi:.15f}"
        print(data_line)
        data['corr_len_ratio_flow'].write(data_line + "\n")

        ######################
        # ent_corr_len_ratio #
        ######################

        data_line = f"{np.log(xi):.15f}\t{SvN:.15f}\t{r:.15f}"
        print(data_line)
        data['ent_corr_len_ratio'].write(data_line + "\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    prog_args, stem_args, leaf_args = fa.parse_input_arguments("ratio_flow")

    my_ratio_flow(prog_args['path'], prog_args['threads'], stem_args['model'], stem_args['chi_max'], leaf_args)
