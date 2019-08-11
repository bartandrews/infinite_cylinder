import numpy as np
import time
import sys
import importlib

import functions as f

parameters_module = "parameters.param_" + str(sys.argv[1])
p = importlib.import_module(parameters_module)


def my_Ly_flow(model, lattice, initial_state, tile_unit, chi_max, chi_max_K, t, U, mu, V, Lx, Ly_min, Ly_max, Ly_samp):

    ent_scal_stem = f.file_name_stem("ent_scal", model, lattice, initial_state, tile_unit, chi_max)
    ent_scal_leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s_%s.dat" % (t, U, mu, V, Lx, Ly_min, Ly_max))
    ent_scal_file = "data/ent_scal/" + ent_scal_stem.replace(" ", "_") + ent_scal_leaf
    open(ent_scal_file, "w")
    ent_scal_data = open(ent_scal_file, "a", buffering=1)

    ##################################################################################################################

    for Ly in np.linspace(Ly_min, Ly_max, Ly_samp, dtype=int):

        (E, psi, M) = f.run_iDMRG_pickle("Ly_flow", model, lattice, initial_state, tile_unit, chi_max,
                                         t, U, mu, V, Lx, Ly, p.use_pickle, p.make_pickle)

        ############
        # ent_scal #
        ############

        print("{Ly:d}\t{SvN:.15f}\t{Sinf:.15f}".format(Ly=Ly, SvN=psi.entanglement_entropy()[0],
                                                       Sinf=psi.entanglement_entropy(n=np.inf)[0]))
        ent_scal_data.write("%i\t%.15f\t%.15f\n" % (Ly, psi.entanglement_entropy()[0],
                                           psi.entanglement_entropy(n=np.inf)[0]))

        #################
        # ent_spec_real #
        #################

        ent_spec_real_stem = f.file_name_stem("ent_spec_real", model, lattice, initial_state, tile_unit, chi_max)
        ent_spec_real_leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V, Lx, Ly))
        ent_spec_real_file = "data/ent_spec_real/" + ent_spec_real_stem.replace(" ", "_") + ent_spec_real_leaf
        open(ent_spec_real_file, "w")
        ent_spec_real_data = open(ent_spec_real_file, "a", buffering=1)

        ##############################################################################################################

        spectrum = psi.entanglement_spectrum(by_charge=True)

        print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=len(spectrum[0][0][0])))

        # spectrum[bond][sector][0][0] --> spectrum[bond][sector][0][n] for different charge entries
        for bond in range(0, Lx*Ly):
            for sector in range(0, len(spectrum[bond])):
                for i in range(0, len(spectrum[bond][sector][1])):
                    print("{charge:d}\t{bond:d}\t{spectrum:.15f}".format(charge=spectrum[bond][sector][0][0],
                                                                         bond=bond,
                                                                         spectrum=spectrum[bond][sector][1][i]))
                    ent_spec_real_data.write("%i\t%i\t%.15f\n" % (spectrum[bond][sector][0][0], bond,
                                                                  spectrum[bond][sector][1][i]))

        ################
        # ent_spec_mom #
        ################

        ent_spec_mom_stem = f.file_name_stem("ent_spec_mom", model, lattice, initial_state, tile_unit, chi_max)
        ent_spec_mom_leaf = ("chi_K_%s_t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat" % (chi_max_K, t, U, mu, V, Lx, Ly))
        ent_spec_mom_file = "data/ent_spec_mom/" + ent_spec_mom_stem.replace(" ", "_") + ent_spec_mom_leaf
        open(ent_spec_mom_file, "w")
        ent_spec_mom_data = open(ent_spec_mom_file, "a", buffering=1)

        ##############################################################################################################

        (Un, W, q, ov, trunc_err) = psi.compute_K(perm=M.lat, trunc_par={'chi_max': chi_max_K}, canonicalize=1.e-6,
                                                  verbose=0)

        if np.abs(np.abs(ov)-1) > 0.1:
            print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))
            print('Warning: State is not invariant under the permutation.')
        else:
            print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))

        print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=q.charges.shape[1]))

        # q.to_qflat()[i][0] --> q.to_qflat()[i][n] for different charge entries
        for i in range(len(W)):
            print("{q:d}\t{K:.15f}\t{epsilon:.15f}".format(q=q.to_qflat()[i][0], K=np.angle(W[i])/np.pi,
                                                           epsilon=-np.log(np.abs(W[i]))))
            ent_spec_mom_data.write("%i\t%.15f\t%.15f\n" % (q.to_qflat()[i][0], np.angle(W[i])/np.pi,
                                                            -np.log(np.abs(W[i]))))


if __name__ == '__main__':

    t0 = time.time()

    my_Ly_flow(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.chi_max_K, p.t, p.U, p.mu, p.V, p.Lx,
               Ly_min=6, Ly_max=8, Ly_samp=2)

    print(time.time() - t0)
