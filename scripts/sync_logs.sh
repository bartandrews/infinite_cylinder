#!/bin/bash
# sync_logs.sh --- sync log_observables from dart to laptop (execute on laptop)
#
# Conditions:
# - extraneous log_observables files are deleted on the laptop so that logs/observables is an exact replica of dart
# - directories are synced but extraneous directories are not deleted from laptop

# rsync options
#
# delete = delete extraneous files from destination
# n = dry run (optional)
# v = verbose
# t = preserve modification time
# z = compress files for transfer
# h = human-readable output
# r = recursive (optional)
# e ssh = transfer files over a secure connection

MODELS="BosHofSqu1 FerHofSqu1"
DIR_PATH=/home/bart/PycharmProjects/infinite_cylinder/logs/observables

# dry run
for DIR in ${MODELS}
do
	# shellcheck disable=SC2029
	if ssh dart "[ -d /home/bart/PycharmProjects/infinite_cylinder/logs/observables/${DIR}/ ]"
	then
		echo
		echo ">>> Dry run from bart@dart:${DIR_PATH}/${DIR}/ to bart@bart:${DIR_PATH}/${DIR}/"
		echo
		rsync -nvtzhre ssh --delete dart:${DIR_PATH}/"${DIR}"/ ${DIR_PATH}/"${DIR}"/
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
		if ssh dart "[ -d /home/bart/PycharmProjects/infinite_cylinder/logs/observables/${DIR}/ ]"
		then
			echo
			echo ">>> Actual run from bart@dart:${DIR_PATH}/${DIR}/ to bart@bart:${DIR_PATH}/${DIR}/"
			echo
			rsync -vtzhre ssh --delete dart:${DIR_PATH}/"${DIR}"/ ${DIR_PATH}/"${DIR}"/
		fi
	done
else
    exit
fi