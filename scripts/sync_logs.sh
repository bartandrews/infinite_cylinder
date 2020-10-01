#!/bin/bash
# sync_logs.sh --- sync log_observables from dirac to laptop (execute on laptop)
#
# Conditions:
# - extraneous log_observables files are deleted on the laptop so that logs/observables is an exact replica of dirac

# rsync options
#
# n = dry run (optional)
# v = verbose
# t = preserve modification time
# z = compress files for transfer
# h = human-readable output
# r = recursive (optional)
# e ssh = transfer files over a secure connection

MODELS="BosHofSqu1 FerHofSqu1"
DIR_PATH_DIRAC=/disk/data11/tfp/BartMadhav/project1/logs/observables
DIR_PATH_BART=/home/bart/PycharmProjects/infinite_cylinder/logs/observables

# dry run
for DIR in ${MODELS}
do
	# shellcheck disable=SC2029
	if ssh bandrews@dirac "[ -d /disk/data11/tfp/BartMadhav/project1/logs/observables/${DIR}/ ]"
	then
		echo
		echo ">>> Dry run from bandrews@dirac:${DIR_PATH_DIRAC}/${DIR}/ to bart@bart:${DIR_PATH_BART}/${DIR}/"
		echo
		rsync -nvtzhre ssh dirac:${DIR_PATH_DIRAC}/"${DIR}"/ ${DIR_PATH_BART}/"${DIR}"/
	fi
done

echo
read -r -p ">>> Are you sure that you want to continue with rsync? [y/N] " response

# actual run
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
    for DIR in ${MODELS}
	do
		# shellcheck disable=SC2029
		if ssh bandrews@dirac "[ -d /disk/data11/tfp/BartMadhav/project1/logs/observables/${DIR}/ ]"
		then
			echo
			echo ">>> Actual run from bandrews@dirac:${DIR_PATH_DIRAC}/${DIR}/ to bart@bart:${DIR_PATH_BART}/${DIR}/"
			echo
			rsync -vtzhre ssh dirac:${DIR_PATH_DIRAC}/"${DIR}"/ ${DIR_PATH_BART}/"${DIR}"/
		fi
	done
else
    exit
fi