# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.file_proc as fp
import functions.dmrg as fd


def my_Ly_flow(model, chi_max, chi_max_K, t1, t2, t2dash, U, mu, V,
               nnvalue, nd_min, nd_max, pvalue, q_min, q_max, nu_samp,
               Lx_MUC, Ly_min, Ly_max, Ly_samp, tag,
               use_pickle, make_pickle):

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U}_mu_{mu}_V_{V}_n_{nnvalue}_{nd_min}_{nd_max}_{nu_samp}_" \
           f"nphi_{pvalue}_{q_min}_{q_max}_{nu_samp}_Lx_MUC_{Lx_MUC}_Ly_{Ly_min}_{Ly_max}_{Ly_samp}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("Ly_flow", model, leaf)

    tools = ["ent_scal", "ent_spec_real", "ent_spec_mom"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf, chi_max_K)

    ####################################################################################################################

    for ndvalue, qvalue in zip(np.linspace(nd_min, nd_max, nu_samp, dtype=int),
                               np.linspace(q_min, q_max, nu_samp, dtype=int)):
        for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

            (data_line, LylB) = fp.print_LylB_headings(model, Ly, ndvalue, nd_min, pvalue, qvalue, q_min, Ly_min)
            if nu_samp != 1 or Ly_samp != 1:
                for tool in ['ent_spec_real', 'ent_spec_mom']:
                    data[tool].write(data_line)

            (E, psi, M) = fd.my_iDMRG_pickle("Ly_flow", model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue,
                                             pvalue, qvalue, Lx_MUC, Ly, use_pickle, make_pickle, run=True)

            ############
            # ent_scal #
            ############

            data_line = "{LylB:.15f}\t{SvN:.15f}\t{Sinf:.15f}"\
                .format(LylB=LylB,
                        SvN=psi.entanglement_entropy()[0],
                        Sinf=psi.entanglement_entropy(n=np.inf)[0])
            print(data_line)
            data['ent_scal'].write(data_line+"\n")

            #################
            # ent_spec_real #
            #################

            spectrum = psi.entanglement_spectrum(by_charge=True)

            print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=len(spectrum[0][0][0])))

            # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
            for bond in range(0, Lx_MUC*Ly):
                for sector in range(0, len(spectrum[bond])):
                    for i in range(0, len(spectrum[bond][sector][1])):
                        data_line = "{charge:d}\t{bond:d}\t{spectrum:.15f}"\
                            .format(charge=spectrum[bond][sector][0][0],
                                    bond=bond,
                                    spectrum=spectrum[bond][sector][1][i])
                        print(data_line)
                        data['ent_spec_real'].write(data_line+"\n")

            ################
            # ent_spec_mom #
            ################

            (Un, W, q, ov, trunc_err) = \
                psi.compute_K(perm=M.lat, trunc_par={'chi_max': chi_max_K}, canonicalize=1.e-6, verbose=0)

            if np.abs(np.abs(ov)-1) > 0.1:
                print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))
                print('Warning: State is not invariant under the permutation.')
            else:
                print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))

            print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=q.charges.shape[1]))

            # q.to_qflat()[i][0] --> q.to_qflat()[i][n] for different charge entries
            for i in range(len(W)):
                data_line = "{q:d}\t{K:.15f}\t{epsilon:.15f}"\
                    .format(q=q.to_qflat()[i][0], K=np.angle(W[i])/np.pi, epsilon=-np.log(np.abs(W[i])))
                print(data_line)
                data['ent_spec_mom'].write(data_line+"\n")

            ##################################
            # corr_func (under construction) #
            ##################################

            if "Orbital" in model:
                op = "Ntot"
            else:
                op = "N"

            # corr_func can be computed beyond psi.L, however then the mps2lat function will not work
            NN = psi.correlation_function(op, op, sites1=range(0, psi.L), sites2=[0])[:, 0]
            print(NN)
            # NN_reshaped = M.lat.mps2lat_values(NN)
            # print(NN_reshaped.shape, M.lat.shape)
            # import pdb; pdb.set_trace()
            # return

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    prc.mkl_set_nthreads(1)

    my_Ly_flow(model="BosonicHofstadter", chi_max=50, chi_max_K=500,
               t1=1, t2=0, t2dash=0, U=0, mu=0, V=0,
               nnvalue=1, nd_min=8, nd_max=8, pvalue=1, q_min=4, q_max=4, nu_samp=1,
               Lx_MUC=1, Ly_min=4, Ly_max=4, Ly_samp=1, tag="",
               use_pickle=False, make_pickle=False)
