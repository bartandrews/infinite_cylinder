# --- python imports
import time
import sys
import numpy as np
# --- TeNPy imports
import tenpy.tools.process as prc
# --- infinite_cylinder imports
import functions.func_proc as fp
import functions.func_dmrg as fd
import functions.func_args as fa


def my_phi_flow(path_flag, threads, model, chi_max, ham_params):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    prc.mkl_set_nthreads(threads)
    t0 = time.time()

    leaf = fp.file_name_leaf("phi_flow", model, ham_params)
    sys.stdout = sys.stderr = fp.Logger("phi_flow", path, model, chi_max, leaf)

    tools = ["overlap", "charge_pump", "ent_spec_flow"]
    data = fp.prepare_output_files(tools, path, model, chi_max, leaf)

    ####################################################################################################################

    engine = fd.my_iDMRG_pickle("phi_flow", path, model, chi_max, ham_params, run=False)

    for phi in np.linspace(ham_params['phi_min'], ham_params['phi_max'], ham_params['phi_samp']):

        if phi == ham_params['phi_min']:
            engine.run()
        else:
            engine.options['mixer'] = False
            del engine.options['chi_list']  # comment out this line for single site DMRG tests
            ham_params.update(phi=phi)
            M = fd.define_iDMRG_model(model, ham_params)
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

    prog_args, stem_args, leaf_args = fa.parse_input_arguments("phi_flow")

    my_phi_flow(prog_args['path'], prog_args['threads'], stem_args['model'], stem_args['chi_max'], leaf_args)
