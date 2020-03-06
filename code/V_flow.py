# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_V_flow(threads, model, chi_max, t1, t2, t2dash, U, mu, V_min, V_max, V_samp,
              nnvalue, nd_min, nd_max, pvalue, q_min, q_max, nu_samp,
              Lx_MUC, Ly_min, Ly_max, Ly_samp, tag,
              use_pickle, make_pickle):

    prc.mkl_set_nthreads(threads)

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U}_mu_{mu}_V_{V_min}_{V_max}_{V_samp}_" \
           f"n_{nnvalue}_{nd_min}_{nd_max}_{nu_samp}_nphi_{pvalue}_{q_min}_{q_max}_{nu_samp}_" \
           f"Lx_MUC_{Lx_MUC}_Ly_{Ly_min}_{Ly_max}_{Ly_samp}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("V_flow", model, leaf)

    tools = ["corr_len_V_flow", "ent_spec_V_flow"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf)

    ##################################################################################################################

    for ndvalue, qvalue in zip(np.linspace(nd_min, nd_max, nu_samp, dtype=int),
                               np.linspace(q_min, q_max, nu_samp, dtype=int)):
        for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

            if nu_samp != 1 or Ly_samp != 1:
                (data_line, _) = fp.print_LylB_headings(model, Ly, ndvalue, nd_min, pvalue, qvalue, q_min, Ly_min)
                for tool in tools:
                    data[tool].write(data_line)

            for V in np.linspace(V_min, V_max, V_samp):
                (E, psi, M) = fd.my_iDMRG_pickle("V_flow", model, chi_max, t1, t2, t2dash, U, mu, V,
                                                 nnvalue, ndvalue, pvalue, qvalue,
                                                 Lx_MUC, Ly, use_pickle, make_pickle, run=True)

                ###################
                # corr_len_V_flow #
                ###################

                xi = psi.correlation_length()

                data_line = f"{V:.15f}\t{xi:.15f}"
                print(data_line)
                data['corr_len_V_flow'].write(data_line+"\n")

                ###################
                # ent_spec_V_flow #
                ###################

                # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
                spectrum = psi.entanglement_spectrum(by_charge=True)

                bond = 0

                for sector in range(0, len(spectrum[bond])):
                    for i in range(0, len(spectrum[bond][sector][1])):
                        data_line = "{charge:d}\t{V:.15f}\t{spectrum:.15f}"\
                            .format(charge=spectrum[bond][sector][0][0],
                                    V=V, spectrum=spectrum[bond][sector][1][i])
                        print(data_line)
                        data['ent_spec_V_flow'].write(data_line+"\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_V_flow(threads=1, model="FerHofHex1Hex5Orbital", chi_max=50,
              t1=1, t2=-0.025, t2dash=0.1, U=0, mu=0, V_min=0, V_max=10, V_samp=100,
              nnvalue=1, nd_min=9, nd_max=9, pvalue=1, q_min=3, q_max=3, nu_samp=1,
              Lx_MUC=1, Ly_min=6, Ly_max=6, Ly_samp=1, tag="",
              use_pickle=False, make_pickle=False)
