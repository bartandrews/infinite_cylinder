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


def my_kappa_flow(path_flag, threads, model, chi_max, ham_params):

    path = "/home/bart/Desktop" if path_flag else ""  # specify the custom path
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("kappa_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("kappa_flow", path, model, chi_max, leaf)

    tools = ["corr_len_kappa_flow", "ent_spec_kappa_flow"]
    data = fp.prepare_output_files(tools, path, model, chi_max, leaf)

    ####################################################################################################################

    t5, t5dash = ham_params['t5'], ham_params['t5dash']  # flow variables

    for kappa in np.linspace(ham_params['kappa_min'], ham_params['kappa_max'], ham_params['kappa_samp']):

        ham_params.update(t5=kappa*t5)
        ham_params.update(t5dash=kappa*t5dash)

        (E, psi, M, _, _) = fd.my_iDMRG_pickle("kappa_flow", path, model, chi_max, ham_params, run=True)

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

    prog_args, stem_args, leaf_args = fa.parse_input_arguments("kappa_flow")

    my_kappa_flow(prog_args['path'], prog_args['threads'], stem_args['model'], stem_args['chi_max'], leaf_args)
