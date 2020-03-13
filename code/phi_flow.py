# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd


def my_phi_flow(threads, model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                nn, nd_min, nd_max, p, q_min, q_max, nu_samp,
                LxMUC, Ly_min, Ly_max, Ly_samp,
                phi_min, phi_max, phi_samp, tag,
                use_pickle, make_pickle):

    prc.mkl_set_nthreads(threads)

    t0 = time.time()

    leaf = f"t1_{t1}_t2_{t2}_t2dash_{t2dash}_U_{U}_mu_{mu}_V_{V}_{Vtype}_{Vrange}_" \
           f"n_{nn}_{nd_min}_{nd_max}_{nu_samp}_nphi_{p}_{q_min}_{q_max}_{nu_samp}_" \
           f"LxMUC_{LxMUC}_Ly_{Ly_min}_{Ly_max}_{Ly_samp}_phi_{phi_min}_{phi_max}_{phi_samp}.dat{tag}"
    sys.stdout = sys.stderr = fp.Logger("phi_flow", model, leaf)

    tools = ["overlap", "charge_pump", "ent_spec_flow"]
    data = fp.prepare_output_files(tools, model, chi_max, leaf)

    ####################################################################################################################

    for nd, q in zip(np.linspace(nd_min, nd_max, nu_samp, dtype=int),
                               np.linspace(q_min, q_max, nu_samp, dtype=int)):
        for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

            if nu_samp != 1 or Ly_samp != 1:
                (data_line, _) = fp.print_LylB_headings(model, Ly, nd, nd_min, p, q, q_min, Ly_min)
                for tool in tools:
                    data[tool].write(data_line)

            engine = fd.my_iDMRG_pickle("phi_flow", model, chi_max, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                                        nn, nd, p, q, LxMUC, Ly,
                                        use_pickle, make_pickle, phi_min, run=False)

            for phi in np.linspace(phi_min, phi_max, phi_samp):

                if phi == phi_min:
                    engine.run()
                else:
                    engine.engine_params['mixer'] = False
                    del engine.engine_params['chi_list']  # comment out this line for single site DMRG tests
                    M = fd.define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, Vtype, Vrange,
                                              nn, nd, p, q, LxMUC, Ly, phi)
                    psi_old = engine.psi
                    engine.init_env(model=M)
                    engine.run()

                    ###########
                    # overlap #
                    ###########

                    abs_ov = abs(psi_old.overlap(engine.psi))
                    data_line = f"{phi:.15f}\t{abs_ov:.15f}"
                    print(data_line)
                    data['overlap'].write(data_line+"\n")

                ###############
                # charge_pump #
                ###############

                QL = engine.psi.average_charge(bond=0)[0]

                data_line = f"{phi:.15f}\t{QL:.15f}"
                print(data_line)
                data['charge_pump'].write(data_line+"\n")

                #################
                # ent_spec_flow #
                #################

                # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
                spectrum = engine.psi.entanglement_spectrum(by_charge=True)

                bond = 0

                for sector in range(0, len(spectrum[bond])):
                    for i in range(0, len(spectrum[bond][sector][1])):
                        data_line = "{charge:d}\t{phi:.15f}\t{spectrum:.15f}"\
                            .format(charge=spectrum[bond][sector][0][0],
                                    phi=phi,
                                    spectrum=spectrum[bond][sector][1][i])
                        print(data_line)
                        data['ent_spec_flow'].write(data_line+"\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == '__main__':

    my_phi_flow(threads=1, model="BosHofSqu1", chi_max=50,
                t1=1, t2=0, t2dash=0, U=0, mu=0, V=0, Vtype='Coulomb', Vrange=0,
                nn=1, nd_min=8, nd_max=8, p=1, q_min=4, q_max=4, nu_samp=1,
                LxMUC=1, Ly_min=4, Ly_max=4, Ly_samp=1, phi_min=0, phi_max=2, phi_samp=21, tag="",
                use_pickle=False, make_pickle=False)
