# --- python imports
import time
import sys
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_ground_state(threads, model, chi_max, use_pickle, **ham_params):

    fp.check_input_params("ground_state", threads, model, chi_max, ham_params)
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("ground_state", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("ground_state", model, leaf)

    ####################################################################################################################

    fd.my_iDMRG_pickle("ground_state", model, chi_max, ham_params, use_pickle, make_pickle=True, run=True)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_ground_state(threads=1, model="BosHofSqu1", chi_max=100, use_pickle=False,
                    t1=1, t2=0, t2dash=0, U=0, mu=0, V=0, Vtype='Coulomb', Vrange=0, n=(1, 8), nphi=(1, 4),
                    LxMUC=1, Ly=4, tag="")
