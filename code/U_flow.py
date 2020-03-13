# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_U_flow(threads, model, chi_max, t1, t2, t2dash, U_min, U_max, U_samp, mu, V, Vtype, Vrange,
              nn, nd_min, nd_max, p, q_min, q_max, nu_samp,
              LxMUC, Ly_min, Ly_max, Ly_samp, tag,
              use_pickle, make_pickle):

    prc.mkl_set_nthreads(threads)

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U_min}_{U_max}_{U_samp}_mu_{mu}_V_{V}_{Vtype}_{Vrange}_" \
           f"n_{nn}_{nd_min}_{nd_max}_{nu_samp}_nphi_{p}_{q_min}_{q_max}_{nu_samp}_" \
           f"LxMUC_{LxMUC}_Ly_{Ly_min}_{Ly_max}_{Ly_samp}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("U_flow", model, leaf)

    tools = ["corr_len_U_flow", "double_occ_U_flow"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf)

    ####################################################################################################################

    for nd, q in zip(np.linspace(nd_min, nd_max, nu_samp, dtype=int),
                               np.linspace(q_min, q_max, nu_samp, dtype=int)):
        for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

            if nu_samp != 1 or Ly_samp != 1:
                (data_line, _) = fp.print_LylB_headings(model, Ly, nd, nd_min, p, q, q_min, Ly_min)
                for tool in tools:
                    data[tool].write(data_line)

            for U in np.linspace(U_min, U_max, U_samp):
                (E, psi, M) = fd.my_iDMRG_pickle("U_flow", model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                                                 nn, nd, p, q,
                                                 LxMUC, Ly, use_pickle, make_pickle, run=True)

                ###################
                # corr_len_U_flow #
                ###################

                xi = psi.correlation_length()

                data_line = f"{U:.15f}\t{xi:.15f}"
                print(data_line)
                data['corr_len_U_flow'].write(data_line+"\n")

                #####################
                # double_occ_U_flow #
                #####################

                term1 = psi.expectation_value('Ntot Ntot')
                term2 = -2*psi.expectation_value('Ntot')

                term1_2 = term1 + term2

                nd_summand = [x + 1 for x in term1_2]
                nd = np.average(nd_summand)

                data_line = f"{U:.15f}\t{nd:.15f}"
                print(data_line)
                data['double_occ_U_flow'].write(data_line+"\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_U_flow(threads=1, model="FerHofHex1Hex5Orbital", chi_max=50,
              t1=1, t2=-0.025, t2dash=0.1, U_min=0, U_max=10, U_samp=100, mu=0, V=0, Vtype='Coulomb', Vrange=1,
              nn=1, nd_min=9, nd_max=9, p=1, q_min=3, q_max=3, nu_samp=1,
              LxMUC=1, Ly_min=6, Ly_max=6, Ly_samp=1, tag="",
              use_pickle=False, make_pickle=False)
