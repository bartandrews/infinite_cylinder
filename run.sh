#!/bin/bash

########################################################################################################################
#
# nohup = no hang-up signal, keeps the job running after you close PyCharm/terminal
# code/... = we need to launch run.sh from the project directory and target a program in the source directory
# >/dev/null 2>&1 & = we send the stdout (>) to null and also the stderr (2>) to the same place (&1).
# (we have our own logging functionality, so we don't need this output)
# We want to run other jobs before this job finishes and so we end with a &
#
########################################################################################################################

#for chi_val in 51 52 53 54
#do
#    nohup python code/phi_flow.py -thr 1 -mod "BosHofSqu1" -chi ${chi_val} -t1 1 -n 1 8 -nphi 1 4 -LxMUC 1 -Ly 4 -phi_min 0 -phi_max 2 -phi_samp 21 >/dev/null 2>&1 &
#done

nohup python code/phi_flow.py -thr 1 -mod "BosHofSqu1" -chi 50 -t1 1 -n 1 8 -nphi 1 4 -LxMUC 1 -Ly 4 -phi_min 0 -phi_max 2 -phi_samp 21 >/dev/null 2>&1 &
