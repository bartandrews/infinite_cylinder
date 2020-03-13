# --- python imports
import time
import sys
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_ground_state(threads, model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                    nn, nd, p, q, LxMUC, Ly, tag):

    prc.mkl_set_nthreads(threads)

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U}_mu_{mu}_V_{V}_{Vtype}_{Vrange}_" \
           f"n_{nn}_{nd}_nphi_{p}_{q}_LxMUC_{LxMUC}_Ly_{Ly}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("ground_state", model, leaf)

    ####################################################################################################################

    fd.my_iDMRG_pickle("ground_state", model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                       nn, nd, p, q, LxMUC, Ly, use_pickle=False, make_pickle=True, run=True)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_ground_state(threads=1, model="BosHofSqu1", chi_max=250,
                    t1=1, t2=0, t2dash=0, U=0, mu=0, V=0, Vtype='Coulomb', Vrange=0,
                    nn=1, nd=8, p=1, q=4,
                    LxMUC=1, Ly=4, tag="")
