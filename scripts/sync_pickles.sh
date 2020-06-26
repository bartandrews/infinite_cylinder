#!/bin/bash
# sync_pickles.sh --- sync pickles from remotes to dart (execute on laptop)
#
# Conditions:
# - username == ubuntu is assumed for all remote hostnames
# - it is assumed that the home directory path is the same on all remote hosts
# - the pickles on remote computers are deleted after they have been successfully copied

# rsync options
#
# remove-source-files = remove the source files after they have been successfully synced
# n = dry run (optional)
# v = verbose
# t = preserve modification time
# z = compress files for transfer
# h = human-readable output
# r = recursive (optional)
# e ssh = transfer files over a secure connection

MODELS="BosHofSqu1 FerHofSqu1"  # the model directories that we would like to sync
PCS="baandr1 baandr2 baandr3"  # remote hostnames that we intend to run on (username == ubuntu is assumed throughout)
UBUNTU_BASE=/home/ubuntu  # the home directory on the remote computers with username == ubuntu
BANDREWS_BASE=/home/fkp/bandrews  # the home directory on the remote computers with username == bandrews
BAANDR_BASE=/home/cluster/baandr  # the home directory on the remote computers with username == baandr
DST_BASE=/home/bart  # the home directory on dart
DIR_PATH=PycharmProjects/infinite_cylinder/pickles/ground_state  # path to the ground_state directory relative to home

# dry run
for DIR in ${MODELS}
do
	for PC in ${PCS}
	do
		# shellcheck disable=SC2029
		if ssh ubuntu@"${PC}" "[ -d /home/ubuntu/PycharmProjects/infinite_cylinder/pickles/ground_state/${DIR}/ ]"
		then
			echo
			echo ">>> Dry run from ubuntu@${PC}:${UBUNTU_BASE}/${DIR_PATH}/${DIR}/ to bart@dart:${DST_BASE}/${DIR_PATH}/${DIR}/"
			echo
			ssh -A dart rsync --remove-source-files -nvtzhre ssh ubuntu@"${PC}":${UBUNTU_BASE}/${DIR_PATH}/"${DIR}"/ ${DST_BASE}/${DIR_PATH}/"${DIR}"/
		fi
	done
	# shellcheck disable=SC2029
	if ssh bandrews@dirac "[ -d /home/fkp/bandrews/PycharmProjects/infinite_cylinder/pickles/ground_state/${DIR}/ ]"
	then
		echo
		echo ">>> Dry run from bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH}/${DIR}/ to bart@dart:${DST_BASE}/${DIR_PATH}/${DIR}/"
		echo
		ssh -A dart rsync --remove-source-files -nvtzhre ssh bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH}/"${DIR}"/ ${DST_BASE}/${DIR_PATH}/"${DIR}"/
	fi
	# shellcheck disable=SC2029
	if ssh baandr@s3it "[ -d /home/cluster/baandr/data/pickles/ground_state/${DIR}/ ]"
	then
		echo
		echo ">>> Dry run from baandr@s3it:${BAANDR_BASE}/data/pickles/ground_state/${DIR}/ to bart@dart:${DST_BASE}/${DIR_PATH}/${DIR}/"
		echo
		ssh -A dart rsync --remove-source-files -nvtzhre ssh baandr@s3it:${BAANDR_BASE}/data/pickles/ground_state/"${DIR}"/ ${DST_BASE}/${DIR_PATH}/"${DIR}"/
	fi
done

echo
read -r -p ">>> Are you sure that you want to continue with rsync? [y/N] " response

# actual run
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
    for DIR in ${MODELS}
	do
		for PC in ${PCS}
		do
			# shellcheck disable=SC2029
			if ssh ubuntu@"${PC}" "[ -d /home/ubuntu/PycharmProjects/infinite_cylinder/pickles/ground_state/${DIR}/ ]"
			then
				echo
				echo ">>> Actual run from ubuntu@${PC}:${UBUNTU_BASE}/${DIR_PATH}/${DIR}/ to bart@dart:${DST_BASE}/${DIR_PATH}/${DIR}/"
				echo
				ssh -A dart rsync --remove-source-files -vtzhre ssh ubuntu@"${PC}":${UBUNTU_BASE}/${DIR_PATH}/"${DIR}"/ ${DST_BASE}/${DIR_PATH}/"${DIR}"/
			fi
		done
		# shellcheck disable=SC2029
		if ssh bandrews@dirac "[ -d /home/fkp/bandrews/PycharmProjects/infinite_cylinder/pickles/ground_state/${DIR}/ ]"
		then
			echo
			echo ">>> Actual run from bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH}/${DIR}/ to bart@dart:${DST_BASE}/${DIR_PATH}/${DIR}/"
			echo
			ssh -A dart rsync --remove-source-files -vtzhre ssh bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH}/"${DIR}"/ ${DST_BASE}/${DIR_PATH}/"${DIR}"/
		fi
		# shellcheck disable=SC2029
		if ssh baandr@s3it "[ -d /home/cluster/baandr/data/pickles/ground_state/${DIR}/ ]"
		then
			echo
			echo ">>> Actual run from baandr@s3it:${BAANDR_BASE}/data/pickles/ground_state/${DIR}/ to bart@dart:${DST_BASE}/${DIR_PATH}/${DIR}/"
			echo
			ssh -A dart rsync --remove-source-files -vtzhre ssh baandr@s3it:${BAANDR_BASE}/data/pickles/ground_state/"${DIR}"/ ${DST_BASE}/${DIR_PATH}/"${DIR}"/
		fi
	done
else
    exit
fi