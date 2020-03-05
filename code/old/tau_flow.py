import numpy as np
import time
import tenpy.tools.process as prc

import functions as f


def my_tau_flow(model, chi_max, t1, t2, t2dash, tau_min, tau_max, tau_samp, U, mu, V, nnvalue, nd_min, nd_max, pvalue, q_min, q_max, nu_samp, Lx, Ly_min, Ly_max, Ly_samp, tag, use_pickle, make_pickle):

    corr_len_tau_flow_stem = f.file_name_stem("corr_len_tau_flow", model, chi_max)
    ent_spec_tau_flow_stem = f.file_name_stem("ent_spec_tau_flow", model, chi_max)
    leaf = ("t1_%s_t2_%s_t2dash_%s_tau_%s_%s_%s_U_%s_mu_%s_V_%s_n_%s_%s_%s_%s_nphi_%s_%s_%s_%s_Lx_%s_Ly_%s_%s_%s.dat%s" % (t1, t2, t2dash, tau_min, tau_max, tau_samp, U, mu, V, nnvalue, nd_min, nd_max, nu_samp, pvalue, q_min, q_max, nu_samp, Lx, Ly_min, Ly_max, Ly_samp, tag))
    corr_len_tau_flow_file = "data/corr_len_tau_flow/" + model + "/" + corr_len_tau_flow_stem.replace(" ", "_") + leaf
    ent_spec_tau_flow_file = "data/ent_spec_tau_flow/" + model + "/" + ent_spec_tau_flow_stem.replace(" ", "_") + leaf
    open(corr_len_tau_flow_file, "w")
    open(ent_spec_tau_flow_file, "w")
    corr_len_tau_flow_data = open(corr_len_tau_flow_file, "a", buffering=1)
    ent_spec_tau_flow_data = open(ent_spec_tau_flow_file, "a", buffering=1)

    ##################################################################################################################

    # engine = f.define_iDMRG_engine_pickle("tau_flow", model, chi_max, t1, t2, t2dash, tau_min*U, mu, tau_min*V, Lx, Ly, use_pickle, make_pickle)

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

            corr_len_tau_flow_data.write(data_line)
            ent_spec_tau_flow_data.write(data_line)

            for tau in np.linspace(tau_min, tau_max, tau_samp):

                (E, psi, M) = f.run_iDMRG_pickle("kappa_flow", model, chi_max, t1, t2, t2dash, tau*U, mu, tau*V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, use_pickle, make_pickle)

                #####################
                # corr_len_tau_flow #
                #####################

                xi = psi.correlation_length()

                data_line = "{tau:.15f}\t{xi:.15f}".format(tau=tau, xi=xi)
                print(data_line)
                corr_len_tau_flow_data.write(data_line+"\n")

                #####################
                # ent_spec_tau_flow #
                #####################

                # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
                spectrum = psi.entanglement_spectrum(by_charge=True)

                bond = 0

                for sector in range(0, len(spectrum[bond])):
                    for i in range(0, len(spectrum[bond][sector][1])):
                        data_line = "{charge:d}\t{tau:.15f}\t{spectrum:.15f}".format(charge=spectrum[bond][sector][0][0], tau=tau, spectrum=spectrum[bond][sector][1][i])
                        print(data_line)
                        ent_spec_tau_flow_data.write(data_line+"\n")


if __name__ == '__main__':

    prc.mkl_set_nthreads(1)

    t0 = time.time()

    my_tau_flow(model="FermionicHex1Hex5Orbital", chi_max=50,
                t1=1, t2=-0.025, t2dash=0.1, tau_min=0, tau_max=1, tau_samp=11, U=100, mu=0, V=10,
                nnvalue=1, nd_min=15, nd_max=15, pvalue=1, q_min=3, q_max=3, nu_samp=1,
                Lx=1, Ly_min=5, Ly_max=5, Ly_samp=1, tag=".polarized2",
                use_pickle=False, make_pickle=False)

    print(time.time() - t0)
