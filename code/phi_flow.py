import numpy as np
import time
import tenpy.tools.process as prc
import sys

import functions as f


class Logger(object):
    def __init__(self, model, leaf):
        self.terminal = sys.stdout or sys.stderr
        self.log = open("data/output/phi_flow/" + model + "/" + "output_phi_flow_" + model + "_" + leaf, 'w', buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


def my_phi_flow(model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, nd_min, nd_max, pvalue, q_min, q_max, nu_samp, Lx_MUC, Ly_min, Ly_max, Ly_samp, phi_min, phi_max, phi_samp, tag, use_pickle, make_pickle):

    overlap_stem = f.file_name_stem("overlap", model, chi_max)
    charge_pump_stem = f.file_name_stem("charge_pump", model, chi_max)
    ent_spec_flow_stem = f.file_name_stem("ent_spec_flow", model, chi_max)
    leaf = ("t1_%s_t2_%s_t2dash_%s_U_%s_mu_%s_V_%s_n_%s_%s_%s_%s_nphi_%s_%s_%s_%s_Lx_MUC_%s_Ly_%s_%s_%s_phi_%s_%s_%s.dat%s" % (t1, t2, t2dash, U, mu, V, nnvalue, nd_min, nd_max, nu_samp, pvalue, q_min, q_max, nu_samp, Lx_MUC, Ly_min, Ly_max, Ly_samp, phi_min, phi_max, phi_samp, tag))
    sys.stdout = sys.stderr = Logger(model, leaf)
    overlap_file = "data/overlap/" + model + "/" + overlap_stem.replace(" ", "_") + leaf
    charge_pump_file = "data/charge_pump/" + model + "/" + charge_pump_stem.replace(" ", "_") + leaf
    ent_spec_flow_file = "data/ent_spec_flow/" + model + "/" + ent_spec_flow_stem.replace(" ", "_") + leaf
    open(overlap_file, "w")
    open(charge_pump_file, "w")
    open(ent_spec_flow_file, "w")
    overlap_data = open(overlap_file, "a", buffering=1)
    charge_pump_data = open(charge_pump_file, "a", buffering=1)
    ent_spec_flow_data = open(ent_spec_flow_file, "a", buffering=1)

    ##################################################################################################################

    for ndvalue, qvalue in zip(np.linspace(nd_min, nd_max, nu_samp, dtype=int), np.linspace(q_min, q_max, nu_samp, dtype=int)):
        for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

            if model in ["BosonicHofstadter", "FermionicHofstadter"]:
                LylB = Ly * np.sqrt(2 * np.pi * (pvalue / qvalue))
            else:
                LylB = Ly * np.sqrt((4 * np.pi * (pvalue / qvalue)) / np.sqrt(3))

            if ndvalue == nd_min and qvalue == q_min and Ly == Ly_min:
                data_line = "LylB={LylB:.15f}\n".format(LylB=LylB)
            else:
                data_line = "\n\nLylB={LylB:.15f}\n".format(LylB=LylB)

            overlap_data.write(data_line)
            charge_pump_data.write(data_line)
            ent_spec_flow_data.write(data_line)

            engine = f.define_iDMRG_engine_pickle("phi_flow", model, chi_max, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, use_pickle, make_pickle, phi_min)

            for phi_ext in np.linspace(phi_min, phi_max, phi_samp):

                if phi_ext == phi_min:
                    engine.run()
                else:
                    engine.engine_params['mixer'] = False
                    del engine.engine_params['chi_list']  # comment out this line for single site DMRG tests
                    M = f.define_iDMRG_model(model, t1, t2, t2dash, U, mu, V, nnvalue, ndvalue, pvalue, qvalue, Lx_MUC, Ly, phi_ext)
                    psi_old = engine.psi
                    engine.init_env(model=M)
                    engine.run()
                    abs_ov = abs(psi_old.overlap(engine.psi))
                    data_line = "{phi_ext:.15f}\t{abs_ov:.15f}".format(phi_ext=phi_ext, abs_ov=abs_ov)
                    print(data_line)
                    overlap_data.write(data_line+"\n")

                ###############
                # charge_pump #
                ###############

                QL = engine.psi.average_charge(bond=0)[0]

                data_line = "{phi_ext:.15f}\t{QL:.15f}".format(phi_ext=phi_ext, QL=QL)
                print(data_line)
                charge_pump_data.write(data_line+"\n")

                #################
                # ent_spec_flow #
                #################

                # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
                spectrum = engine.psi.entanglement_spectrum(by_charge=True)

                bond = 0

                for sector in range(0, len(spectrum[bond])):
                    for i in range(0, len(spectrum[bond][sector][1])):
                        data_line = "{charge:d}\t{phi_ext:.15f}\t{spectrum:.15f}".format(charge=spectrum[bond][sector][0][0], phi_ext=phi_ext, spectrum=spectrum[bond][sector][1][i])
                        print(data_line)
                        ent_spec_flow_data.write(data_line+"\n")


if __name__ == '__main__':

    prc.mkl_set_nthreads(1)

    t0 = time.time()

    # my_phi_flow(model="BosonicHofstadter", chi_max=50,
    #             t1=1, t2=0, t2dash=0, U=0, mu=0, V=0,
    #             nnvalue=1, nd_min=8, nd_max=8, pvalue=1, q_min=4, q_max=4, nu_samp=1,
    #             Lx_MUC=1, Ly_min=4, Ly_max=4, Ly_samp=1, phi_min=0, phi_max=2, phi_samp=21, tag="",
    #             use_pickle=False, make_pickle=False)

    my_phi_flow(model="FermionicHex1Hex5Orbital", chi_max=50,
                t1=1, t2=-0.0125, t2dash=0.05, U=100, mu=0, V=10,
                nnvalue=1, nd_min=9, nd_max=9, pvalue=1, q_min=3, q_max=3, nu_samp=1,
                Lx_MUC=1, Ly_min=6, Ly_max=6, Ly_samp=1, phi_min=0, phi_max=3, phi_samp=31, tag="",
                use_pickle=False, make_pickle=False)

    print(time.time() - t0)
