# --- python imports
import time
import sys
import pickle
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_obser as fo
import functions.func_args as fa


def my_observables(pickle_file, threads, chiK_max):

    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    (model, chi_max, leaf, LxMUC, Ly, extra_dof_flag) = fp.process_pickle_file_name(pickle_file)
    sys.stdout = sys.stderr = fp.Logger("observables", model, chi_max, leaf)

    # Here, you need to enter the tools that you are interested in studying.
    tools = ["ent_spec_real", "ent_spec_mom", "corr_func"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf, chiK_max)

    ####################################################################################################################

    with open(pickle_file, 'rb') as file1:
        [E, psi, M, _, _] = pickle.load(file1)

    fo.scalar_observables(E, psi)
    # fo.nonscalar_observables(tools, data, psi, M, chiK_max, LxMUC, Ly, extra_dof_flag, print_data=False)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    (file, prog_args, obser_args) = fa.parse_observables_input_arguments()

    my_observables(file, prog_args['threads'], obser_args['chiK_max'])
