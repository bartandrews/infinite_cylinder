# --- python imports
import time
import sys
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_ground_state(threads, model, chi_max, t1, t2, t2dash, U, mu, V,
                    nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, tag):

    prc.mkl_set_nthreads(threads)

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U}_mu_{mu}_V_{V}_n_{nnvalue}_{ndvalue}_" \
           f"nphi_{pvalue}_{qvalue}_Lx_MUC_{Lx_MUC}_Ly_{Ly}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("ground_state", model, leaf)

    ####################################################################################################################

    fd.my_iDMRG_pickle("ground_state", model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue,
                        pvalue, qvalue, Lx_MUC, Ly, use_pickle=False, make_pickle=True, run=True)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_ground_state(threads=1, model="BosHofSqu1", chi_max=250,
                    t1=1, t2=0, t2dash=0, U=0, mu=0, V=0,
                    nnvalue=1, ndvalue=8, pvalue=1, qvalue=4,
                    Lx_MUC=1, Ly=4, tag="")
