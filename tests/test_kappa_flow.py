from kappa_flow import my_kappa_flow


def test_my_kappa_flow():

    my_kappa_flow(threads=1, model="FermionicHex1Hex5Orbital", chi_max=150,
                  t1=1, t2=-0.025, t2dash=0.1, kappa_min=0, kappa_max=1, kappa_samp=11, U=100, mu=0, V=10,
                  nnvalue=1, nd_min=9, nd_max=9, pvalue=1, q_min=3, q_max=3, nu_samp=1,
                  Lx_MUC=1, Ly_min=6, Ly_max=6, Ly_samp=1, tag="",
                  use_pickle=False, make_pickle=False)

    return
