#!/bin/bash

########################################################################################################################
# --- REFERENCE ---
#
# nohup = no hang-up signal, keeps the job running after you close PyCharm/terminal
# nice = run the code with a niceness level of 10 (you can specifiy more precisely with -n flag)
# code/... = we need to launch run.sh from the project directory and target a program in the source directory
# >/dev/null 2>&1 = we send the stdout (>) to null and also the stderr (2>) to the same place (&1).
# (we have our own logging functionality, so we don't need this output)
# If we want the following jobs to start before the current job finishes, we end the line with a &
#
########################################################################################################################
# --- EXAMPLES ---
#
#for chi_val in 51 52 53 54
#do
#    nohup nice python code/ground_state.py -thr 1 -mod "FerHofSqu1" -chi ${chi_val} -t1 1 -V 10 -Vtype "Coulomb" -Vrange 1 -n 1 9 -nphi 1 3 -LxMUC 1 -Ly 6 >/dev/null 2>&1 &
#done
#
# nohup nice python code/ground_state.py -thr 1 -mod "BosHofSqu1" -chi 10 -t1 1 -n 1 8 -nphi 1 4 -LxMUC 1 -Ly 4 >/dev/null 2>&1
# nohup nice python code/observables.py -thr 1 -chiK 250 pickles/ground_state/BosHofSqu1/E_psi_M_BosHofSqu1_chi_10_t1_1_n_1_8_nphi_1_4_LxMUC_1_Ly_4.pkl >/dev/null 2>&1 &
#
########################################################################################################################

nohup PYTHONPATH=~/TeNPy python code/ground_state.py -thr 1 -mod "FerHofSqu1" -chi 500 -t1 1 -V 10 -Vtype "Coulomb" -Vrange 1 -n 1 9 -nphi 1 3 -LxMUC 1 -Ly 9 #>/dev/null 2>&1 &

#nohup python code/observables.py -thr 1 -chiK 500 pickles/ground_state/FerHofSqu1/E_psi_M_FerHofSqu1_chi_500_t1_1_V_10_Vtype_Coulomb_Vrange_1_n_1_9_nphi_1_3_LxMUC_1_Ly_9.pkl >/dev/null 2>&1 &
