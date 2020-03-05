from V_flow import my_V_flow


def test_my_V_flow():

    my_V_flow(threads=1, model="FermionicHex1Hex5Orbital", chi_max=50,
              t1=1, t2=-0.025, t2dash=0.1, U=0, mu=0, V_min=0, V_max=10, V_samp=100,
              nnvalue=1, nd_min=9, nd_max=9, pvalue=1, q_min=3, q_max=3, nu_samp=1,
              Lx_MUC=1, Ly_min=6, Ly_max=6, Ly_samp=1, tag="",
              use_pickle=False, make_pickle=False)

    return
