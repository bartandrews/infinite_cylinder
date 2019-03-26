import numpy as np
import time

import functions as f
import parameters as p


def my_ent_spec_mom(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly, charge_sectors):

    stem = f.file_name_stem("ent_spec_mom", model, lattice, initial_state, tile_unit, chi_max, charge_sectors)
    leaf = ("t_%s_U_%s_mu_%s_V_%s_Lx_%s_Ly_%s.dat" % (t, U, mu, V, Lx, Ly))
    dat_file = stem + leaf
    open(dat_file, "w")
    data = open(dat_file, "a", buffering=1)

    (E, psi, M) = f.run_iDMRG(model, lattice, initial_state, tile_unit, chi_max, t, U, mu, V, Lx, Ly)

    (Un, W, q, ov, trunc_err) = psi.compute_K(perm=M.lat, canonicalize=1.e-6, verbose=0)

    if np.abs(np.abs(ov)-1) > 0.1:
        print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))
        print('Warning: State is not invariant under the permutation.')
    else:
        print("|ov|={ov_abs:.15f}".format(ov_abs=np.abs(ov)))

    if charge_sectors:

        print("We select charge entry 1 out of qnumber={qnumber:d}.".format(qnumber=q.charges.shape[1]))

        # q.to_qflat()[i][0] --> q.to_qflat()[i][n] for different charge entries
        for i in range(len(W)):
            print("{q:d}\t{K:.15f}\t{epsilon:.15f}".format(q=q.to_qflat()[i][0], K=np.angle(W[i])/np.pi,
                                                           epsilon=-np.log(np.abs(W[i]))))
            data.write("%i\t%.15f\t%.15f\n" % (q.to_qflat()[i][0], np.angle(W[i])/np.pi, -np.log(np.abs(W[i]))))

    else:

        for i in range(len(W)):
            print("{K:.15f}\t{epsilon:.15f}".format(K=np.angle(W[i]) / np.pi, epsilon=-np.log(np.abs(W[i]))))
            data.write("%.15f\t%.15f\n" % (np.angle(W[i]) / np.pi, -np.log(np.abs(W[i]))))


if __name__ == '__main__':

    t0 = time.time()

    my_ent_spec_mom(p.model, p.lattice, p.initial_state, p.tile_unit, p.chi_max, p.t, p.U, p.mu, p.V, p.Lx, p.Ly,
                    charge_sectors=True)

    print(time.time() - t0)
