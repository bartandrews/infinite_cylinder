# --- python imports
import time
import sys
import pickle
import gzip
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_obser as fo
import functions.func_args as fa


def my_observables(pickle_file, path_flag, threads, scalar, chiK_max):

    path = "/home/bart/Desktop" if path_flag else ""  # specify the custom path
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    (model, chi_max, leaf, extra_dof_flag) = fp.process_pickle_file_name(pickle_file)
    sys.stdout = sys.stderr = fp.Logger("observables", path, model, chi_max, leaf)

    if not scalar:
        # Here, you need to enter the tools that you are interested in studying.
        tools = ["ent_spec_real", "ent_spec_mom", "density", "corr_func"]
        data = fp.prepare_output_files(tools, path, model, chi_max, leaf, chiK_max)

    ####################################################################################################################

    with (gzip.open if fp.is_gz_file(pickle_file) else open)(pickle_file, 'rb') as file1:
        state_data = pickle.load(file1)

    E, psi, M = state_data['E'], state_data['psi'], state_data['M']

    fo.scalar_observables(E, psi)
    if not scalar:
        fo.nonscalar_observables(tools, data, psi, M, chiK_max, extra_dof_flag, print_data=False)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    (file, prog_args, obser_args) = fa.parse_observables_input_arguments()

    my_observables(file, prog_args['path'], prog_args['threads'], prog_args['scalar'], obser_args['chiK_max'])
