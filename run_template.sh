#!/bin/bash

# --- Definitions needed for PyCharm remote SSH external tools --- #####################################################
#
#TENPY_DIR=~/TeNPy
#SOURCE_DIR=~/PycharmProjects/infinite_cylinder/code
#
#export PYTHONPATH=${TENPY_DIR}:${SOURCE_DIR}
#
#PYTHON_EXE=/usr/local/anaconda3/bin/python3  # for [bart, dart]
##PYTHON_EXE=~/anaconda3/bin/python3  # for [baandr1, baandr2, baandr3]
#
## ampersand (&) still does not work...
#
# --- Reference --- ####################################################################################################
#
# nohup = no hang-up signal, keeps the job running after you close PyCharm/terminal
# nice = run the code with a niceness level of 10 (you can specifiy more precisely with -n flag)
# parallel = keep a certain number of jobs running at any one time (default = number of cores)
# code/... = we need to launch run.sh from the project directory and target a program in the source directory
# >/dev/null 2>&1 = we send the stdout (>) to null and also the stderr (2>) to the same place (&1).
# (we have our own logging functionality, so we don't need this output)
# > nohup.out = this overwrites the nohup.out file on each run, and only stores the output of the preceeding program
# e.g. parallel (writes which jobs are yet to finish after 1 kill signal, etc.)
# If we want the following jobs to start before the current job finishes, we end the line with a &
# = sign after an optional argument accounts for possible negative sign ambiguities
#
# --- Example --- ######################################################################################################

# original hopping parameters
t1_full=1
t5_full=-0.025
t5dash_full=0.1

########
# bart #
########

if [ "$HOSTNAME" == "bart" ]
then
    shopt -s nullglob  # account for the case of no .pkl files in the path
    FILES="pickles/ground_state/FerHofHex1Hex5Orbital/*.pkl"
    for file in $FILES
    do
        python code/observables.py -thr 4 -chiK 500 "${file}"
    done
fi

########
# dart #
########

if [ "$HOSTNAME" == "dart" ]
then
    # charge pumping for the nu=2/5 state
    charge_pump() {
    for chi_val in {50..150..50}
    do
        for i in {0..10}
        do
            t5_val=$(bc -l <<<"$t5_full * ($i / 10)" | awk '{printf "%g\n", $0}')
            t5dash_val=$(bc -l <<<"$t5dash_full * ($i / 10)" | awk '{printf "%g\n", $0}')
            echo python code/phi_flow.py -thr 1 -mod "FerHofHex1Hex5Orbital" -chi ${chi_val} -t1 $t1_full -t5="${t5_val}" -t5dash="${t5dash_val}" -U 100 -V 10 -Vtype "Coulomb" -Vrange 1 -n 1 15 -nphi 1 3 -LxMUC 1 -Ly 5 -phi_min 0 -phi_max 5 -phi_samp 51
        done
    done
    }
    export -f charge_pump
    charge_pump | nohup parallel -j 8 > nohup.out &
fi

###########
# baandr1 #
###########

if [ "$HOSTNAME" == "baandr1" ]
then
    # charge pumping for the nu=1/3 state
    charge_pump() {
    for chi_val in {50..150..50}
    do
        for i in {0..10}
        do
            t5_val=$(bc -l <<<"$t5_full * ($i / 10)" | awk '{printf "%g\n", $0}')
            t5dash_val=$(bc -l <<<"$t5dash_full * ($i / 10)" | awk '{printf "%g\n", $0}')
            echo python code/phi_flow.py -thr 1 -mod "FerHofHex1Hex5Orbital" -chi ${chi_val} -t1 $t1_full -t5="${t5_val}" -t5dash="${t5dash_val}" -U 100 -V 10 -Vtype "Coulomb" -Vrange 1 -n 1 18 -nphi 1 3 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 3 -phi_samp 31
        done
    done
    }
    export -f charge_pump
    charge_pump | nohup parallel > nohup.out &
fi

###########
# baandr2 #
###########

if [ "$HOSTNAME" == "baandr2" ]
then
    # entanglement scaling for the nu=1/3 state
    ent_scal() {
    for chi_val in {300..500..100}
    do
        for Ly_val in 6 9
        do
            for q_val in 3 4 5
            do
                for i in 0 3 4
                do
                    t5_val=$(bc -l <<<"$t5_full * ($i / 10)" | awk '{printf "%g\n", $0}')
                    t5dash_val=$(bc -l <<<"$t5dash_full * ($i / 10)" | awk '{printf "%g\n", $0}')
                    nd_val_sixth=$(bc -l <<<"6 * $q_val" | awk '{printf "%g\n", $0}')
                    echo python code/ground_state.py -thr 1 -mod "FerHofHex1Hex5Orbital" -chi ${chi_val} -t1 $t1_full -t5="${t5_val}" -t5dash="${t5dash_val}" -U 100 -V 10 -Vtype "Coulomb" -Vrange 1 -n 1 "${nd_val_sixth}" -nphi 1 ${q_val} -LxMUC 1 -Ly ${Ly_val}
                done
            done
        done
    done
    }
    export -f ent_scal
    ent_scal | nohup parallel > nohup.out &
fi

###########
# baandr3 #
###########

if [ "$HOSTNAME" == "baandr3" ]
then
    # kappa flow for the nu=1/3 and nu=2/5 states
    kappa_flow() {
    for chi_val in {50..150..50}
    do
        echo python code/kappa_flow.py -thr 1 -mod "FerHofHex1Hex5Orbital" -chi ${chi_val} -t1 $t1_full -t5=$t5_full -t5dash=$t5dash_full -kappa_min 0 -kappa_max 1 -kappa_samp 11 -U 100 -V 10 -Vtype "Coulomb" -Vrange 1 -n 1 18 -nphi 1 3 -LxMUC 1 -Ly 6
        echo python code/kappa_flow.py -thr 1 -mod "FerHofHex1Hex5Orbital" -chi ${chi_val} -t1 $t1_full -t5=$t5_full -t5dash=$t5dash_full -kappa_min 0 -kappa_max 1 -kappa_samp 11 -U 100 -V 10 -Vtype "Coulomb" -Vrange 1 -n 1 15 -nphi 1 3 -LxMUC 1 -Ly 5
    done
    }
    export -f kappa_flow
    kappa_flow | nohup parallel > nohup.out &
fi
