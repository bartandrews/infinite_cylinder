#!/bin/bash
# run.sh --- run programs on remotes (execute on remotes)
#
# Conditions:
# - this script needs to be executed from the project root directory, since relative paths are used for program execution

########
# dart #
########

if [ "$HOSTNAME" == "dart" ]
then
	runs() {  # CI MWEs (Hofstadter)
	for chi_val in 50 100 200 300 400 500
	do
		echo python code/phi_flow.py -thr 4 -mod "FerHofSqu1" -chi ${chi_val} -t1 1 -n 1 14 -nphi 3 14 -LxMUC 1 -Ly 9 -phi_min 0 -phi_max 1 -phi_samp 11
		echo python code/phi_flow.py -thr 4 -mod "FerHofSqu1" -chi ${chi_val} -t1 1 -n 1 19 -nphi 4 19 -LxMUC 1 -Ly 12 -phi_min 0 -phi_max 1 -phi_samp 11
	done
    }
    export -f runs
    runs | nohup parallel -j 2 > nohup.out &
fi

###########
# baandr1 #
###########

if [ "$HOSTNAME" == "baandr1" ]
then
    runs() {  # FCI MWEs (Haldane)
	for chi_val in 50 100 200 300 400 500
	do
		echo python code/phi_flow.py -thr 2 -mod "FerHalSquC1" -chi ${chi_val} -t1 1 -V 2 -Vtype Coulomb -Vrange 1 -n 1 6 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 3 -phi_samp 31
		echo python code/phi_flow.py -thr 2 -mod "FerHalSquC1" -chi ${chi_val} -t1 1 -V 2 -Vtype Coulomb -Vrange 1 -n 1 6 -LxMUC 1 -Ly 9 -phi_min 0 -phi_max 3 -phi_samp 31
	done
    }
    export -f runs
    runs | nohup parallel -j 4 > nohup.out &
fi

###########
# baandr2 #
###########

if [ "$HOSTNAME" == "baandr2" ]
then
	runs() {  # FCI Haldane
	for chi_val in 50 100 200 300 400 500
	do
		echo python code/phi_flow.py -thr 2 -mod "FerHalSquC1" -chi ${chi_val} -t1 1 -V 2 -Vtype Coulomb -Vrange 1 -n 1 6 -LxMUC 1 -Ly 12 -phi_min 0 -phi_max 3 -phi_samp 31
	done
    }
    export -f runs
    runs | nohup parallel -j 4 > nohup.out &
fi

###########
# baandr3 #
###########

if [ "$HOSTNAME" == "baandr3" ]
then
	runs() {  # CI MWEs (Hofstadter) Ly=6
	for chi_val in 50 100 200 300 400 500
	do
		echo python code/phi_flow.py -thr 1 -mod "FerHofSqu1" -chi ${chi_val} -t1 1 -n 4 9 -nphi 4 9 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 1 -phi_samp 11
		echo python code/phi_flow.py -thr 1 -mod "FerHofSqu1" -chi ${chi_val} -t1 1 -n 5 14 -nphi 5 14 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 1 -phi_samp 11
		echo python code/phi_flow.py -thr 1 -mod "FerHofSqu1" -chi ${chi_val} -t1 1 -n 5 19 -nphi 5 19 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 1 -phi_samp 11
		echo python code/phi_flow.py -thr 1 -mod "FerHofSqu1" -chi ${chi_val} -t1 1 -n 5 24 -nphi 5 24 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 1 -phi_samp 11

	done
    }
    export -f runs
    runs | nohup parallel -j 8 > nohup.out &
fi

############
# baandr4 #
############

if [ "$HOSTNAME" == "baandr4" ]
then
	runs() {
	for chi_val in 100 200 400 600
    do
		echo python code/phi_flow.py -thr 8 -mod "BosHalC1" -chi ${chi_val} -t1 1 -n 1 1 -LxMUC 1 -Ly 4 -phi_min 0 -phi_max 1 -phi_samp 11
		echo python code/phi_flow.py -thr 8 -mod "BosHalC1" -chi ${chi_val} -t1 1 -n 1 2 -LxMUC 1 -Ly 4 -phi_min 0 -phi_max 2 -phi_samp 21
	done
    }
    export -f runs
    runs | nohup parallel -j 2 > nohup.out &
fi

############
# baandr5 #
############

if [ "$HOSTNAME" == "baandr5" ]
then
	runs() {
	for chi_val in 100 200 400 600
    do
		echo python code/phi_flow.py -thr 16 -mod "FerHalC1" -chi ${chi_val} -t1 1 -n 1 3 -V 1 -Vtype "Coulomb" -Vrange 1 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 3 -phi_samp 31
	done
    }
    export -f runs
    runs | nohup parallel -j 1 > nohup.out &
fi

############
# baandr6 #
############

if [ "$HOSTNAME" == "baandr6" ]
then
	runs() {
	for chi_val in 100 200 400 600
    do
		echo python code/phi_flow.py -thr 8 -mod "BosHalC3" -chi ${chi_val} -t1 1 -n 1 1 -LxMUC 1 -Ly 4 -phi_min 0 -phi_max 1 -phi_samp 11
		echo python code/phi_flow.py -thr 8 -mod "BosHalC3" -chi ${chi_val} -t1 1 -n 1 2 -LxMUC 1 -Ly 4 -phi_min 0 -phi_max 2 -phi_samp 21
	done
    }
    export -f runs
    runs | nohup parallel -j 2 > nohup.out &
fi

############
# baandr7 #
############

if [ "$HOSTNAME" == "baandr7" ]
then
	runs() {
	for chi_val in 100 200 400 600
    do
		echo python code/phi_flow.py -thr 8 -mod "FerHalC3" -chi ${chi_val} -t1 1 -n 1 1 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 1 -phi_samp 11
		echo python code/phi_flow.py -thr 8 -mod "FerHalC3" -chi ${chi_val} -t1 1 -n 1 3 -V 1 -Vtype "Coulomb" -Vrange 1 -LxMUC 1 -Ly 6 -phi_min 0 -phi_max 3 -phi_samp 31
	done
    }
    export -f runs
    runs | nohup parallel -j 2 > nohup.out &
fi

#########
# dirac #
#########

if [ "$HOSTNAME" == "dirac" ]
then
	runs() {
#		# correcting fermions at nu=1/3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 2 7 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 2 7 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 2 7 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 2 7 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1450 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 2 7 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1050 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 2 7 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 250 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 500 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 250 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 500 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1150 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1150 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1150 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1150 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 24 -nphi 3 8 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 30 -nphi 3 10 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 30 -nphi 3 10 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 30 -nphi 3 10 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 30 -nphi 3 10 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 30 -nphi 3 10 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 30 -nphi 3 10 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 30 -nphi 3 10 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 30 -nphi 3 10 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 5
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 5
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 11
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 11
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 3 11 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 250 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 2
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 500 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 2
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 5
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 5
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 11
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 11
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 39 -nphi 3 13 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 5
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 3 14 -LxMUC 1 -Ly 11
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 5
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 11
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 3 17 -LxMUC 1 -Ly 11
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 57 -nphi 3 19 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 57 -nphi 3 19 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 57 -nphi 3 19 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 57 -nphi 3 19 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 4 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 33 -nphi 4 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 4 15 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 4 15 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 4 15 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 4 15 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 4 17 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 4 17 -LxMUC 1 -Ly 3
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 4 17 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 4 17 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 4 17 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 4 17 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 6 17 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 51 -nphi 6 17 -LxMUC 1 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 1 27 -nphi 2 9 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1000 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 1 33 -nphi 2 11 -LxMUC 1 -Ly 9
		# correcting fermions at nu=2/5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 35 -nphi 2 7 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 35 -nphi 2 7 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 35 -nphi 2 7 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 35 -nphi 2 7 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 35 -nphi 2 7 -LxMUC 1 -Ly 15
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 850 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 35 -nphi 2 7 -LxMUC 1 -Ly 15
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 45 -nphi 2 9 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 45 -nphi 2 9 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 45 -nphi 2 9 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 45 -nphi 2 9 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 45 -nphi 2 9 -LxMUC 1 -Ly 15
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 55 -nphi 2 11 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 55 -nphi 2 11 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 850 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 55 -nphi 2 13 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 900 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 55 -nphi 2 13 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 20 -nphi 3 8 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 20 -nphi 3 8 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 25 -nphi 3 10 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 25 -nphi 3 10 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 850 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 55 -nphi 3 11 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 900 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 55 -nphi 3 11 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 65 -nphi 5 13 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 65 -nphi 5 13 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1550 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 65 -nphi 5 13 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1600 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 65 -nphi 5 13 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 5 14 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 5 14 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 5 14 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 5 14 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 5 14 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 5 14 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 200 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 5 14 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 5 14 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 40 -nphi 5 16 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 40 -nphi 5 16 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 40 -nphi 5 16 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 40 -nphi 5 16 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 40 -nphi 5 16 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 40 -nphi 5 16 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 40 -nphi 5 16 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 40 -nphi 5 16 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 85 -nphi 5 17 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 85 -nphi 5 17 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 85 -nphi 5 17 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 85 -nphi 5 17 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 85 -nphi 5 17 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 5 18 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 5 18 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 5 18 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 5 18 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 5 18 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 5 18 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 45 -nphi 5 18 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 95 -nphi 5 19 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 95 -nphi 5 19 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 95 -nphi 5 19 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 95 -nphi 5 19 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 95 -nphi 5 19 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 95 -nphi 5 19 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 105 -nphi 5 21 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 105 -nphi 5 21 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 105 -nphi 5 21 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 105 -nphi 5 21 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 2 105 -nphi 5 21 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 2 55 -nphi 2 11 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 2 55 -nphi 2 11 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 2 65 -nphi 2 13 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 2 65 -nphi 2 13 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2950 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 2 55 -nphi 2 11 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 3000 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 2 55 -nphi 2 11 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 2 65 -nphi 2 13 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 2 65 -nphi 2 13 -LxMUC 1 -Ly 10
		# correcting fermions at nu=3/7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 2 9 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 2 9 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 77 -nphi 2 11 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 77 -nphi 2 11 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 91 -nphi 2 13 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 91 -nphi 2 13 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 91 -nphi 2 13 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 91 -nphi 2 13 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 2 15 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 35 -nphi 2 15 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 2 19 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 2 19 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 56 -nphi 3 8 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 56 -nphi 3 8 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 3250 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 70 -nphi 3 10 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 3300 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 70 -nphi 3 10 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2400 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 77 -nphi 3 11 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2450 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 77 -nphi 3 11 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 91 -nphi 3 13 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 91 -nphi 3 13 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 98 -nphi 3 14 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 98 -nphi 3 14 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 77 -nphi 4 11 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 77 -nphi 4 11 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 91 -nphi 4 13 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 5 18 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 5 18 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 5 19 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 5 19 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 119 -nphi 6 17 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 119 -nphi 6 17 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 6 19 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 6 19 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 4
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 4
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 42 -nphi 7 18 -LxMUC 1 -Ly 10
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 7 19 -LxMUC 1 -Ly 4
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 7 19 -LxMUC 1 -Ly 4
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 7 19 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 7 19 -LxMUC 1 -Ly 5
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 7 19 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 133 -nphi 7 19 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 140 -nphi 7 20 -LxMUC 1 -Ly 4
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 140 -nphi 7 20 -LxMUC 1 -Ly 4
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 140 -nphi 7 20 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 800 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 140 -nphi 7 20 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 750 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 140 -nphi 7 20 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 3 77 -nphi 2 11 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 3 77 -nphi 2 11 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 3 91 -nphi 2 13 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 2 -n 3 91 -nphi 2 13 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 3 77 -nphi 2 11 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 3 77 -nphi 2 11 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 3 91 -nphi 2 13 -LxMUC 1 -Ly 14
		echo python code/ground_state.py -thr 1 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 3 -n 3 91 -nphi 2 13 -LxMUC 1 -Ly 14
	}
    export -f runs
    runs | nohup parallel -j 40 > nohup.out &
fi

########
# weyl #
########

if [ "$HOSTNAME" == "weyl" ]
then
	runs() {
#		# correcting bosons at nu=1/2
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 14 -nphi 2 7 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 18 -nphi 2 9 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 22 -nphi 2 11 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 22 -nphi 2 11 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 22 -nphi 2 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 22 -nphi 2 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 2 13 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 16 -nphi 3 8 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 16 -nphi 3 8 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 1 16 -nphi 3 8 -LxMUC 2 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 1 16 -nphi 3 8 -LxMUC 2 -Ly 7
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 1 16 -nphi 3 8 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 1 16 -nphi 3 8 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 1 16 -nphi 3 8 -LxMUC 2 -Ly 9
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 20 -nphi 3 10 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 20 -nphi 3 10 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 20 -nphi 3 10 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 20 -nphi 3 10 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 1 22 -nphi 3 11 -LxMUC 2 -Ly 11
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 3 13 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 3 13 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 3 13 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 3 13 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 3 13 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 3 13 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 28 -nphi 3 14 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 28 -nphi 3 14 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 28 -nphi 3 14 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 28 -nphi 3 14 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 28 -nphi 3 14 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 28 -nphi 3 14 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 28 -nphi 3 14 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 28 -nphi 3 14 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 34 -nphi 3 17 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 34 -nphi 3 17 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 34 -nphi 3 17 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 34 -nphi 3 17 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 34 -nphi 3 17 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 34 -nphi 3 17 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 34 -nphi 3 17 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 34 -nphi 3 17 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 8
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 26 -nphi 4 13 -LxMUC 1 -Ly 12
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 30 -nphi 4 15 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 30 -nphi 4 15 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 30 -nphi 4 15 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 30 -nphi 4 15 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 34 -nphi 4 17 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 34 -nphi 4 17 -LxMUC 1 -Ly 4
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 550 -t1 1 -n 1 34 -nphi 4 17 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 600 -t1 1 -n 1 34 -nphi 4 17 -LxMUC 1 -Ly 6
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 1 28 -nphi 5 14 -LxMUC 1 -Ly 10
#		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 1 38 -nphi 6 19 -LxMUC 1 -Ly 8
		# correcting bosons at nu=2/3
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1200 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 3 -Ly 7
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1200 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 3 -Ly 7
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 1 -Ly 12
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 1 -Ly 12
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 1 -Ly 15
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 21 -nphi 2 7 -LxMUC 1 -Ly 15
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 27 -nphi 2 9 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 27 -nphi 2 9 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 27 -nphi 2 9 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 27 -nphi 2 9 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 27 -nphi 2 9 -LxMUC 1 -Ly 12
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 33 -nphi 2 11 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 33 -nphi 2 11 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 33 -nphi 2 11 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 33 -nphi 2 11 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 1 12 -nphi 3 8 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1200 -t1 1 -n 1 12 -nphi 3 8 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 1 12 -nphi 3 8 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1200 -t1 1 -n 1 12 -nphi 3 8 -LxMUC 1 -Ly 7
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 1 12 -nphi 3 8 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1200 -t1 1 -n 1 12 -nphi 3 8 -LxMUC 1 -Ly 8
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 1 12 -nphi 3 8 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 33 -nphi 3 11 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 33 -nphi 3 11 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 33 -nphi 3 11 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 33 -nphi 3 11 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 33 -nphi 3 11 -LxMUC 1 -Ly 15
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 39 -nphi 3 13 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 39 -nphi 3 13 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 39 -nphi 3 13 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 39 -nphi 3 13 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 1 21 -nphi 3 14 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 1 21 -nphi 3 14 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 1 21 -nphi 3 14 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 1 21 -nphi 3 14 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 2 33 -nphi 4 11 -LxMUC 3 -Ly 7
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 2 33 -nphi 4 11 -LxMUC 3 -Ly 8
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 1150 -t1 1 -n 2 39 -nphi 4 13 -LxMUC 3 -Ly 7
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 45 -nphi 4 15 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 45 -nphi 4 15 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 45 -nphi 4 15 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 45 -nphi 4 15 -LxMUC 1 -Ly 9
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 45 -nphi 4 15 -LxMUC 1 -Ly 15
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 51 -nphi 4 17 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 800 -t1 1 -n 2 51 -nphi 4 17 -LxMUC 1 -Ly 6
		echo python code/ground_state.py -thr 1 -mod BosHofSqu1 -chi 750 -t1 1 -n 2 51 -nphi 4 17 -LxMUC 1 -Ly 9
    }
    export -f runs
    runs | nohup parallel -j 36 > nohup.out &
fi

#######
# new #
#######

if [ "$HOSTNAME" == "new" ]
then
	runs() {
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 49 -nphi 1 7 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 49 -nphi 1 7 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 6 91 -nphi 2 13 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 6 91 -nphi 2 13 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 14 -nphi 1 6 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 14 -nphi 1 6 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 6 77 -nphi 2 11 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 6 77 -nphi 2 11 -LxMUC 1 -Ly 14
	# configurations output
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 9 70 -nphi 3 10 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 9 70 -nphi 3 10 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 12 91 -nphi 4 13 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 12 91 -nphi 4 13 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 13 96 -nphi 6 19 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 13 96 -nphi 6 19 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 9 56 -nphi 3 8 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 9 56 -nphi 3 8 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 12 77 -nphi 4 11 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 12 77 -nphi 4 11 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 13 86 -nphi 6 17 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 13 86 -nphi 6 17 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 15 91 -nphi 5 13 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 15 91 -nphi 5 13 -LxMUC 1 -Ly 7
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 20 -nphi 7 20 -LxMUC 1 -Ly 8
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 20 -nphi 7 20 -LxMUC 1 -Ly 8
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 28 -nphi 1 12 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 28 -nphi 1 12 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 70 -nphi 1 10 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 70 -nphi 1 10 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 19 -nphi 7 19 -LxMUC 1 -Ly 8
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 19 -nphi 7 19 -LxMUC 1 -Ly 8
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 1 9 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 21 -nphi 1 9 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 56 -nphi 1 8 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 56 -nphi 1 8 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 20 -nphi 7 20 -LxMUC 1 -Ly 9
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 20 -nphi 7 20 -LxMUC 1 -Ly 9
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 14 -nphi 1 6 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 14 -nphi 1 6 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 19 -nphi 7 19 -LxMUC 1 -Ly 9
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 19 -nphi 7 19 -LxMUC 1 -Ly 9
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 35 -nphi 1 5 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 3 35 -nphi 1 5 -LxMUC 1 -Ly 14
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 1950 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 6 -nphi 7 18 -LxMUC 1 -Ly 9
	echo python code/ground_state.py -thr 16 -mod FerHofSqu1 -chi 2000 -t1 1 -V 10 -Vtype Coulomb -Vrange 1 -n 1 6 -nphi 7 18 -LxMUC 1 -Ly 9
    }
    export -f runs
    runs | nohup parallel -j 1 > nohup.out &
fi