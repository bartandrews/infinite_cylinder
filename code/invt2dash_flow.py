import numpy as np
import time
import tenpy.tools.process as prc

import functions as f


def my_invt2dash_flow(model, chi_max, t1, t2, invt2dash_min, invt2dash_max, invt2dash_samp, U, mu, V, nnvalue, nd_min, nd_max, pvalue, q_min, q_max, nu_samp, Lx, Ly_min, Ly_max, Ly_samp, tag, use_pickle, make_pickle):

    corr_len_invt2dash_flow_stem = f.file_name_stem("corr_len_invt2dash_flow", model, chi_max)
    ent_spec_invt2dash_flow_stem = f.file_name_stem("ent_spec_invt2dash_flow", model, chi_max)
    leaf = ("t1_%s_t2_%s_invt2dash_%s_%s_%s_U_%s_mu_%s_V_%s_n_%s_%s_%s_%s_nphi_%s_%s_%s_%s_Lx_%s_Ly_%s_%s_%s.dat%s" % (t1, t2, invt2dash_min, invt2dash_max, invt2dash_samp, U, mu, V, nnvalue, nd_min, nd_max, nu_samp, pvalue, q_min, q_max, nu_samp, Lx, Ly_min, Ly_max, Ly_samp, tag))
    corr_len_invt2dash_flow_file = "data/corr_len_invt2dash_flow/" + model + "/" + corr_len_invt2dash_flow_stem.replace(" ", "_") + leaf
    ent_spec_invt2dash_flow_file = "data/ent_spec_invt2dash_flow/" + model + "/" + ent_spec_invt2dash_flow_stem.replace(" ", "_") + leaf
    open(corr_len_invt2dash_flow_file, "w")
    open(ent_spec_invt2dash_flow_file, "w")
    corr_len_invt2dash_flow_data = open(corr_len_invt2dash_flow_file, "a", buffering=1)
    ent_spec_invt2dash_flow_data = open(ent_spec_invt2dash_flow_file, "a", buffering=1)

    ##################################################################################################################

    # engine = f.define_iDMRG_engine_pickle("invt2dash_flow", model, chi_max, t1, t2, 1/invt2dash_min, U, mu, V, Lx, Ly, use_pickle, make_pickle)

    for ndvalue, qvalue in zip(np.linspace(nd_min, nd_max, nu_samp, dtype=int), np.linspace(q_min, q_max, nu_samp, dtype=int)):
        for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

            if model in ["BosonicHofstadter", "FermionicHofstadter"]:
                LylB = Ly * np.sqrt(2*np.pi*(pvalue/qvalue))
            else:
                LylB = Ly * np.sqrt((4*np.pi*(pvalue/qvalue))/np.sqrt(3))

            if ndvalue == nd_min and qvalue == q_min and Ly == Ly_min:
                data_line = "LylB={LylB:.15f}\n".format(LylB=LylB)
            else:
                data_line = "\n\nLylB={LylB:.15f}\n".format(LylB=LylB)

            corr_len_invt2dash_flow_data.write(data_line)
            ent_spec_invt2dash_flow_data.write(data_line)

            for invt2dash in np.linspace(invt2dash_min, invt2dash_max, invt2dash_samp):

                (E, psi, M) = f.run_iDMRG_pickle("invt2dash_flow", model, chi_max, t1, t2, 1/invt2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, use_pickle, make_pickle)

                ###########################
                # corr_len_invt2dash_flow #
                ###########################

                xi = psi.correlation_length()

                data_line = "{invt2dash:.15f}\t{xi:.15f}".format(invt2dash=invt2dash, xi=xi)
                print(data_line)
                corr_len_invt2dash_flow_data.write(data_line+"\n")

                ###########################
                # ent_spec_invt2dash_flow #
                ###########################

                # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
                spectrum = psi.entanglement_spectrum(by_charge=True)

                bond = 0

                for sector in range(0, len(spectrum[bond])):
                    for i in range(0, len(spectrum[bond][sector][1])):
                        data_line = "{charge:d}\t{invt2dash:.15f}\t{spectrum:.15f}".format(charge=spectrum[bond][sector][0][0], invt2dash=invt2dash, spectrum=spectrum[bond][sector][1][i])
                        print(data_line)
                        ent_spec_invt2dash_flow_data.write(data_line+"\n")


if __name__ == '__main__':

    prc.mkl_set_nthreads(1)

    t0 = time.time()

    my_invt2dash_flow(model="BosonicHex1Hex5Orbital", chi_max=50,
                      t1=1, t2=-0.025, invt2dash_min=10, invt2dash_max=30, invt2dash_samp=101, U=0, mu=0, V=0,
                      nnvalue=1, nd_min=8, nd_max=8, pvalue=1, q_min=4, q_max=4, nu_samp=1,
                      Lx=1, Ly_min=6, Ly_max=6, Ly_samp=1, tag="",
                      use_pickle=False, make_pickle=False)

    print(time.time() - t0)
