# --- python imports
import time
import sys
import pickle
import gzip
import os
import h5py
# --- TeNPy imports
import tenpy.tools.process as prc
from tenpy.tools import hdf5_io
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_obser as fo
import functions.func_args as fa


def my_observables(pickle_file_path, path_flag, threads, scalar, chiK_max):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    pickle_file = os.path.split(pickle_file_path)[1]
    (model, chi_max, leaf, extra_dof_flag) = fp.process_pickle(pickle_file)
    sys.stdout = sys.stderr = fp.Logger("observables", path, model, chi_max, leaf)

    if not scalar:
        # Here, you need to enter the tools that you are interested in studying.
        tools = ["ent_spec_real", "ent_spec_mom", "density", "corr_func"]
        data = fp.prepare_output_files(tools, path, model, chi_max, leaf, chiK_max)

    ####################################################################################################################

    if ".h5" in pickle_file:
        with h5py.File(pickle_file_path, 'r') as file1:
            state_data = hdf5_io.load_from_hdf5(file1)
    elif ".pkl" in pickle_file:  # backward compatibility
        with (gzip.open if fp.is_gz_file(pickle_file_path) else open)(pickle_file_path, 'rb') as file1:
            if "E_psi_M" in pickle_file:
                [E, psi, M, _, _] = pickle.load(file1)
            else:
                state_data = pickle.load(file1)
    else:
        raise ValueError("Unknown file extension.")

    if "E_psi_M" not in pickle_file or ".h5" in pickle_file:
        E, psi, M = state_data['E'], state_data['psi'], state_data['M']

    fo.scalar_observables(E, psi)
    if not scalar:
        fo.nonscalar_observables(tools, data, psi, M, chiK_max, extra_dof_flag, print_data=False)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    (file_path, prog_args, obser_args) = fa.parse_observables_input_arguments()

    my_observables(file_path, prog_args['path'], prog_args['threads'], prog_args['scalar'], obser_args['chiK_max'])
