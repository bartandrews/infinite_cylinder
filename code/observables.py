# --- python imports
import time
import sys
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd
import functions.func_obser as fo


def my_observables(threads, model, chi_max, **ham_params):

    fp.check_input_params("observables", threads, model, chi_max, ham_params)
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("observables", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("observables", model, leaf)

    # Here, you need to enter the tools that you are interested in studying.
    tools = ["ent_spec_real", "ent_spec_mom", "corr_func"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf, ham_params['chiK_max'])

    ####################################################################################################################

    # Here, you need to enter the flow that you are targeting to analyze.
    (E, psi, M) = fd.my_iDMRG_pickle("observables", model, chi_max, ham_params,
                                     use_pickle=True, make_pickle=False, run=True)

    fo.scalar_observables(E, psi)
    fo.nonscalar_observables(tools, data, psi, M, ham_params['chiK_max'], ham_params['LxMUC'], ham_params['Ly'],
                             extra_dof=False, print_data=False)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_observables(threads=1, model="BosHofSqu1", chi_max=10, chiK_max=500,
                   t1=1, t2=0, t2dash=0, U=0, mu=0, V=0, Vtype='Coulomb', Vrange=0, n=(1, 8), nphi=(1, 4),
                   LxMUC=1, Ly=4, tag="")
