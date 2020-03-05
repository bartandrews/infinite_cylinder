# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.file_proc as fp
import functions.dmrg as fd


def my_kappa_flow(model, chi_max, t1, t2, t2dash, kappa_min, kappa_max, kappa_samp, U, mu, V,
                  nnvalue, nd_min, nd_max, pvalue, q_min, q_max, nu_samp,
                  Lx_MUC, Ly_min, Ly_max, Ly_samp, tag,
                  use_pickle, make_pickle):

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_kappa_{kappa_min}_{kappa_max}_{kappa_samp}_U_{U}_mu_{mu}_V_{V}_" \
           f"n_{nnvalue}_{nd_min}_{nd_max}_{nu_samp}_nphi_{pvalue}_{q_min}_{q_max}_{nu_samp}_" \
           f"Lx_MUC_{Lx_MUC}_Ly_{Ly_min}_{Ly_max}_{Ly_samp}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("kappa_flow", model, leaf)

    tools = ["corr_len", "ent_spec"]

    stem, file, data = [dict()] * 3
    for tool in tools:
        stem.update({tool: fp.file_name_stem(tool, model, chi_max)})
        file.update({tool: f"data/{tool}/{model}/" + stem[tool].replace(" ", "_") + leaf})
        open(file[tool], "w")
        data[tool] = open(file[tool], "a", buffering=1)

    ####################################################################################################################

    # In case you ever want to update the engine on each step rather than recalculating
    # engine = fd.my_iDMRG_pickle("kappa_flow", model, chi_max, t1, kappa_min*t2, kappa_min*t2dash, U, mu, V,
    #                             nnvalue, nd_min, pvalue, q_min,
    #                             Lx_MUC, Ly_min, use_pickle, make_pickle, run=False)

    for ndvalue, qvalue in zip(np.linspace(nd_min, nd_max, nu_samp, dtype=int),
                               np.linspace(q_min, q_max, nu_samp, dtype=int)):
        for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

            if nu_samp != 1 or Ly_samp != 1:
                (data_line, _) = fp.print_LylB_headings(model, Ly, ndvalue, nd_min, pvalue, qvalue, q_min, Ly_min)
                for tool in tools:
                    data[tool].write(data_line)

            for kappa in np.linspace(kappa_min, kappa_max, kappa_samp):

                (E, psi, M) = fd.my_iDMRG_pickle("kappa_flow", model, chi_max, t1, kappa*t2, kappa*t2dash, U, mu, V,
                                                 nnvalue, ndvalue, pvalue, qvalue,
                                                 Lx_MUC, Ly, use_pickle, make_pickle, run=True)

                #######################
                # corr_len_kappa_flow #
                #######################

                xi = psi.correlation_length()

                data_line = f"{kappa:.15f}\t{xi:.15f}"
                print(data_line)
                data['corr_len_kappa_flow'].write(data_line+"\n")

                #######################
                # ent_spec_kappa_flow #
                #######################

                # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
                spectrum = psi.entanglement_spectrum(by_charge=True)

                bond = 0

                for sector in range(0, len(spectrum[bond])):
                    for i in range(0, len(spectrum[bond][sector][1])):
                        data_line = "{charge:d}\t{kappa:.15f}\t{spectrum:.15f}"\
                            .format(charge=spectrum[bond][sector][0][0],
                                    kappa=kappa, spectrum=spectrum[bond][sector][1][i])
                        print(data_line)
                        data['ent_spec_kappa_flow'].write(data_line+"\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    prc.mkl_set_nthreads(1)

    my_kappa_flow(model="FermionicHex1Hex5Orbital", chi_max=150,
                  t1=1, t2=-0.025, t2dash=0.1, kappa_min=0, kappa_max=1, kappa_samp=11, U=100, mu=0, V=10,
                  nnvalue=1, nd_min=9, nd_max=9, pvalue=1, q_min=3, q_max=3, nu_samp=1,
                  Lx_MUC=1, Ly_min=6, Ly_max=6, Ly_samp=1, tag=".polarized2",
                  use_pickle=False, make_pickle=False)
