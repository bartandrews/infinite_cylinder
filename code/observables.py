# --- python imports
import time
import sys
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd
import functions.func_obser as fo


def my_observables(threads, model, chi_max, chi_max_K, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                   nnvalue, ndvalue, pvalue, qvalue, LxMUC, Ly, tag):

    prc.mkl_set_nthreads(threads)

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U}_mu_{mu}_V_{V}_{Vtype}_{Vrange}_" \
           f"n_{nnvalue}_{ndvalue}_nphi_{pvalue}_{qvalue}_LxMUC_{LxMUC}_Ly_{Ly}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("observables", model, leaf)

    # Here, you need to enter the tools that you are interested in studying.
    tools = ["ent_spec_real", "ent_spec_mom", "corr_func"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf, chi_max_K)

    ####################################################################################################################

    # Here, you need to enter the flow that you are targeting to analyze.
    (E, psi, M) = fd.my_iDMRG_pickle("ground_state", model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                                     nnvalue, ndvalue, pvalue, qvalue, LxMUC, Ly,
                                     use_pickle=True, make_pickle=False, run=True)

    fo.scalar_observables(E, psi)
    fo.nonscalar_observables(tools, data, psi, M, chi_max_K, LxMUC, Ly, extra_dof=False, print_data=True)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_observables(threads=1, model="BosHofSqu1", chi_max=500, chi_max_K=500,
                   t1=1, t2=0, t2dash=0, U=0, mu=0, V=0, Vtype='Coulomb', Vrange=0,
                   nnvalue=1, ndvalue=8, pvalue=1, qvalue=4,
                   LxMUC=1, Ly=4, tag="")
