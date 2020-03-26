#!/bin/bash

#nohup python code/ground_state.py -thr 1 -mod "FerHofSqu1" -chi 500 -t1 1 -V 10 -Vtype "Coulomb" -Vrange 1 -n 1 9 -nphi 1 3 -LxMUC 1 -Ly 9 #>/dev/null 2>&1 &
#nohup python code/phi_flow.py -thr 1 -mod "BosHofSqu1" -chi 51 -t1 1 -n 1 8 -nphi 1 4 -LxMUC 1 -Ly 4 -phi_min 0 -phi_max 2 -phi_samp 21 #>/dev/null 2>&1 &

for chi_val in 51 52 53 54
do
    nohup python code/phi_flow.py -thr 1 -mod "BosHofSqu1" -chi ${chi_val} -t1 1 -n 1 8 -nphi 1 4 -LxMUC 1 -Ly 4 -phi_min 0 -phi_max 2 -phi_samp 21 >/dev/null 2>&1 &
done

#nohup python code/observables.py -thr 1 -chiK 500 pickles/ground_state/FerHofSqu1/E_psi_M_FerHofSqu1_chi_500_t1_1_V_10_Vtype_Coulomb_Vrange_1_n_1_9_nphi_1_3_LxMUC_1_Ly_9.pkl >/dev/null 2>&1 &
