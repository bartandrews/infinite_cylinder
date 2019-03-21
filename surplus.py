# function used to check if Hamiltonian is Hermitian from tests/tes_model.py (insert at the beginning)


import tenpy.linalg.np_conserved as npc
import numpy.testing as npt


def check_hermitian(H):
    """Check if `H` is a hermitian MPO."""
    if not H.finite:
        # include once over the boundary: double the unit cell
        # a general MPO might have terms going over multiple unit cells, but we ignore that...
        Ws = H._W * 2
    else:
        Ws = H._W
    #check trace(H.H) = trace(H.H^dagger)
    W = Ws[0].take_slice([H.get_IdL(0)], ['wL'])

    trHH = npc.tensordot(W, W.replace_label('wR', 'wR*'), axes=[['p', 'p*'], ['p*', 'p']])
    trHHd = npc.tensordot(W, W.conj(), axes=[['p', 'p*'], ['p*', 'p']])
    for W in Ws[1:]:
        trHH = npc.tensordot(trHH, W, axes=['wR', 'wL'])
        trHHd = npc.tensordot(trHHd, W, axes=['wR', 'wL'])
        trHH = npc.tensordot(
            trHH, W.replace_label('wR', 'wR*'), axes=[['wR*', 'p', 'p*'], ['wL', 'p*', 'p']])
        trHHd = npc.tensordot(trHHd, W.conj(), axes=[['wR*', 'p', 'p*'], ['wL*', 'p*', 'p']])
    i = H.get_IdR(H.L - 1)
    trHH = trHH[i, i]
    trHHd = trHHd[i, i]
    print("check_hermitian: ", trHH, trHHd)
    npt.assert_array_almost_equal_nulp(trHH, trHHd, H.L * 20)


# insert after run_iDMRG
check_hermitian(M.H_MPO)


# print the couplings of the model (insert after M is defined)


import matplotlib.pyplot as plt
ax = plt.gca()
M.lat.plot_sites(ax)
M.plot_coupling_terms(ax)
ax.set_aspect("equal")
ax.legend()
plt.show()


# print bond dimension of each bond in the MPO (insert after M is defined)

print(M.H_MPO.chi)


# alternate ways for computing the left charge (use instead of QL_bar and QL)

QL_bar_m = [engine.psi.average_charge(bond=m) for m in range(engine.psi.L)]
QL = (QL_bar_m[0] - np.mean(QL_bar_m, axis=0))[0]
print(phi_ext, QL)