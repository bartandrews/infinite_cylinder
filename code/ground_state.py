# --- python imports
import time
import sys
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd
import functions.func_args as fa


def my_ground_state(threads, model, chi_max, ham_params, use_pickle):

    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("ground_state", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("ground_state", model, chi_max, leaf)

    ####################################################################################################################

    fd.my_iDMRG_pickle("ground_state", model, chi_max, ham_params, use_pickle, make_pickle=True, run=True)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    prog_args, stem_args, leaf_args = fa.parse_input_arguments("ground_state")

    my_ground_state(prog_args['threads'], stem_args['model'], stem_args['chi_max'], leaf_args, prog_args['use_pickle'])
