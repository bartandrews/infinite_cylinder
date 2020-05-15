#!/bin/bash

# This script will run the scalar observables for all new pickles

########
# bart #
########

if [ "$HOSTNAME" == "bart" ]
then
	DIR_PATH="/home/bart/PycharmProjects/infinite_cylinder"
	LOG_PATH="logs/observables"
	PKL_PATH="pickles/ground_state"
	MODELS="BosHofSqu1 FerHofSqu1"

	for MODEL in ${MODELS}
	do
		LOG_DIR="${DIR_PATH}/${LOG_PATH}/${MODEL}"
		PKL_DIR="${DIR_PATH}/${PKL_PATH}/${MODEL}"
		PKL_FILES="${PKL_DIR}/*"

		for PKL_FILE in ${PKL_FILES}
		do
			if test -f "$PKL_FILE"; then
				CONV_FILE_0="$(basename "$PKL_FILE")"  # strip the file from the path
				CONV_FILE_1="${CONV_FILE_0/E_psi_M/log_observables}"  # replace E_psi_M with log_observables (for backward compatibility)
				CONV_FILE_2="${CONV_FILE_1/state/log_observables}"  # replace state with log_observables
				CONV_FILE_3="${CONV_FILE_2/.pkl/.dat}"  # replace .pkl with .dat
				if ! [ -e "$LOG_DIR"/"$CONV_FILE_3" ]; then  # if the corresponding log file is not present
					echo python code/observables.py -s "$PKL_DIR"/"$(basename "$PKL_FILE")"
				fi
			fi
		done
	done
fi

########
# dart #
########

if [ "$HOSTNAME" == "dart" ]
then
	DIR_PATH="/home/bart/PycharmProjects/infinite_cylinder"
	LOG_PATH="logs/observables"
	PKL_PATH="pickles/ground_state"
	MODELS="BosHofSqu1 FerHofSqu1"

	observables() {
	for MODEL in ${MODELS}
	do
		LOG_DIR="${DIR_PATH}/${LOG_PATH}/${MODEL}"
		PKL_DIR="${DIR_PATH}/${PKL_PATH}/${MODEL}"
		PKL_FILES="${PKL_DIR}/*"

		for PKL_FILE in ${PKL_FILES}
		do
			if test -f "$PKL_FILE"; then
				CONV_FILE_0="$(basename "$PKL_FILE")"  # strip the file from the path
				CONV_FILE_1="${CONV_FILE_0/E_psi_M/log_observables}"  # replace E_psi_M with log_observables (for backward compatibility)
				CONV_FILE_2="${CONV_FILE_1/state/log_observables}"  # replace state with log_observables
				CONV_FILE_3="${CONV_FILE_2/.pkl/.dat}"  # replace .pkl with .dat
				# echo $CONV_FILE_3
				if ! [ -e "$LOG_DIR"/"$CONV_FILE_3" ]; then  # if the corresponding log file is not present
					echo python code/observables.py -s "$PKL_DIR"/"$(basename "$PKL_FILE")"
				fi
			fi
		done
	done
	}
    export -f observables
    observables | nohup parallel -j 8 > nohup.out &
fi