import numpy as np
import time
import tenpy.tools.process as prc

import functions as f


def my_Ly_flow(model, chi_max, chi_max_K, t1, t2, t2dash, U, mu, V, nnvalue, nd_min, nd_max, pvalue, q_min, q_max, nu_samp, Lx, Ly_min, Ly_max, Ly_samp, tag, use_pickle, make_pickle):

    ent_scal_stem = f.file_name_stem("ent_scal", model, chi_max)
    ent_spec_real_stem = f.file_name_stem("ent_spec_real", model, chi_max)
    ent_spec_mom_stem = f.file_name_stem("ent_spec_mom", model, chi_max)
    leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_n_%s_%s_%s_%s_nphi_%s_%s_%s_%s_Lx_%s_Ly_%s_%s_%s.dat%s" % (t1, t2, t2dash, U, mu, V, nnvalue, nd_min, nd_max, nu_samp, pvalue, q_min, q_max, nu_samp, Lx, Ly_min, Ly_max, Ly_samp, tag))
    ent_scal_file = "data/ent_scal/" + model + "/" + ent_scal_stem.replace(" ", "_") + leaf
    ent_spec_real_file = "data/ent_spec_real/" + model + "/" + ent_spec_real_stem.replace(" ", "_") + leaf
    ent_spec_mom_file = "data/ent_spec_mom/" + model + "/" + ent_spec_mom_stem.replace(" ", "_") + ("chi_K_%s_" % chi_max_K) + leaf
    open(ent_scal_file, "w")
    open(ent_spec_real_file, "w")
    open(ent_spec_mom_file, "w")
    ent_scal_data = open(ent_scal_file, "a", buffering=1)
    ent_spec_real_data = open(ent_spec_real_file, "a", buffering=1)
    ent_spec_mom_data = open(ent_spec_mom_file, "a", buffering=1)

    ####################################################################################################################

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

            ent_spec_real_data.write(data_line)
            ent_spec_mom_data.write(data_line)

            (E, psi, M) = f.run_iDMRG_pickle("Ly_flow", model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx, Ly, use_pickle, make_pickle)

            ############
            # ent_scal #
            ############

            data_line = "{LylB:.15f}\t{SvN:.15f}\t{Sinf:.15f}".format(LylB=LylB, SvN=psi.entanglement_entropy()[0], Sinf=psi.entanglement_entropy(n=np.inf)[0])
            print(data_line)
            ent_scal_data.write(data_line+"\n")

            #################
            # ent_spec_real #
            #################

            spectrum = psi.entanglement_spectrum(by_charge=True)

            print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=len(spectrum[0][0][0])))

            # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
            for bond in range(0, Lx*Ly):
                for sector in range(0, len(spectrum[bond])):
                    for i in range(0, len(spectrum[bond][sector][1])):
                        data_line = "{charge:d}\t{bond:d}\t{spectrum:.15f}".format(charge=spectrum[bond][sector][0][0], bond=bond, spectrum=spectrum[bond][sector][1][i])
                        print(data_line)
                        ent_spec_real_data.write(data_line+"\n")

            ################
            # ent_spec_mom #
            ################

            (Un, W, q, ov, trunc_err) = psi.compute_K(perm=M.lat, trunc_par={'chi_max': chi_max_K}, canonicalize=1.e-6, verbose=0)

            if np.abs(np.abs(ov)-1) > 0.1:
                print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))
                print('Warning: State is not invariant under the permutation.')
            else:
                print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))

            print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=q.charges.shape[1]))

            # q.to_qflat()[i][0] --> q.to_qflat()[i][n] for different charge entries
            for i in range(len(W)):
                data_line = "{q:d}\t{K:.15f}\t{epsilon:.15f}".format(q=q.to_qflat()[i][0], K=np.angle(W[i])/np.pi, epsilon=-np.log(np.abs(W[i])))
                print(data_line)
                ent_spec_mom_data.write(data_line+"\n")


if __name__ == '__main__':

    prc.mkl_set_nthreads(1)

    t0 = time.time()

    my_Ly_flow(model="FermionicHex1Hex5", chi_max=400, chi_max_K=500,
               t1=1, t2=-0.005, t2dash=0.02, U=100, mu=0, V=10,
               nnvalue=1, nd_min=15, nd_max=15, pvalue=1, q_min=3, q_max=3, nu_samp=1,
               Lx=1, Ly_min=10, Ly_max=10, Ly_samp=1, tag=".polarized2",
               use_pickle=False, make_pickle=False)

    # my_Ly_flow(model="FermionicHex1Hex5Orbital", chi_max=400, chi_max_K=500,
    #            t1=1, t2=-0.01, t2dash=0.04, U=100, mu=0, V=10,
    #            nnvalue=1, nd_min=9, nd_max=9, pvalue=1, q_min=3, q_max=3, nu_samp=1,
    #            Lx=1, Ly_min=9, Ly_max=9, Ly_samp=1, tag="",
    #            use_pickle=False, make_pickle=False)

    print(time.time() - t0)
