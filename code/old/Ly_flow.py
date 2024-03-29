# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd
import functions.func_obser as fo


def my_Ly_flow(threads, model, chi_max, chiK_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
               nn, nd_min, nd_max, p, q_min, q_max, nu_samp,
               LxMUC, Ly_min, Ly_max, Ly_samp, tag,
               use_pickle, make_pickle):

    prc.mkl_set_nthreads(threads)

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U}_mu_{mu}_V_{V}_{Vtype}_{Vrange}_" \
           f"n_{nn}_{nd_min}_{nd_max}_{nu_samp}_nphi_{p}_{q_min}_{q_max}_{nu_samp}_" \
           f"LxMUC_{LxMUC}_Ly_{Ly_min}_{Ly_max}_{Ly_samp}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("Ly_flow", model, leaf)

    tools = ["ent_scal", "ent_spec_real", "ent_spec_mom"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf, chiK_max)

    ####################################################################################################################

    for nd, q in zip(np.linspace(nd_min, nd_max, nu_samp, dtype=int),
                               np.linspace(q_min, q_max, nu_samp, dtype=int)):
        for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

            (data_line, LylB) = fp.print_LylB_headings(model, Ly, nd, nd_min, p, q, q_min, Ly_min)
            if nu_samp != 1 or Ly_samp != 1:
                for tool in ['ent_spec_real', 'ent_spec_mom']:
                    data[tool].write(data_line)

            (E, psi, M) = fd.my_iDMRG_pickle("Ly_flow", model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                                             nn, nd, p, q,
                                             LxMUC, Ly, use_pickle, make_pickle, run=True)

            ############
            # ent_scal #
            ############

            data_line = "{LylB:.15f}\t{SvN:.15f}\t{Sinf:.15f}"\
                .format(LylB=LylB,
                        SvN=psi.entanglement_entropy()[0],
                        Sinf=psi.entanglement_entropy(n=np.inf)[0])
            print(data_line)
            data['ent_scal'].write(data_line+"\n")

            ###############
            # observables #
            ###############

            fo.nonscalar_observables(tools, data, model, psi, M, chiK_max, LxMUC, Ly, print_data=True)

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_Ly_flow(threads=1, model="BosHofSqu1", chi_max=50, chiK_max=500,
               t1=1, t2=0, t2dash=0, U=0, mu=0, V=0, Vtype='Coulomb', Vrange=0,
               nn=1, nd_min=8, nd_max=8, p=1, q_min=4, q_max=4, nu_samp=1,
               LxMUC=1, Ly_min=4, Ly_max=4, Ly_samp=1, tag="",
               use_pickle=False, make_pickle=False)
