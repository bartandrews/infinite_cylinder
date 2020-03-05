from Ly_flow import my_Ly_flow


def test_my_Ly_flow():

    my_Ly_flow(threads=1, model="BosonicHofstadter", chi_max=50, chi_max_K=500,
               t1=1, t2=0, t2dash=0, U=0, mu=0, V=0,
               nnvalue=1, nd_min=8, nd_max=8, pvalue=1, q_min=4, q_max=4, nu_samp=1,
               Lx_MUC=1, Ly_min=4, Ly_max=4, Ly_samp=1, tag="",
               use_pickle=False, make_pickle=False)

    return
