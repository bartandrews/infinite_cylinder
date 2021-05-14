#!/bin/bash
# sync_Vflow.sh --- sync V_flow from various to laptop (execute on laptop)

# rsync options
#
# n = dry run (optional)
# v = verbose
# t = preserve modification time
# z = compress files for transfer
# h = human-readable output
# r = recursive (optional)
# e ssh = transfer files over a secure connection

TOOLS="ent_spec_V_flow ent_V_flow corr_len_V_flow ent_corr_len"
PCS="baandr1 baandr2 baandr3"
BANDREWS_BASE=/home/fkp/bandrews  # the home directory on the remote computers with username == bandrews
BART_BASE=/home/bart  # the home directory on the local computer with username == bart
UBUNTU_BASE=/home/ubuntu  # the home directory on the remote computers with username == ubuntu
DIR_PATH_HOME=PycharmProjects/infinite_cylinder/data  # path to the ground_state directory relative to home

# dry run
for TOOL in ${TOOLS}
do
	# shellcheck disable=SC2029
	if ssh bandrews@dirac "[ -d /home/fkp/bandrews/PycharmProjects/infinite_cylinder/data/${TOOL}/FerHofSqu1/ ]"
	then
		echo
		echo ">>> Dry run from bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/"
		echo
		rsync -nvtzhre ssh bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/
	fi
	# shellcheck disable=SC2029
#	if ssh bart@dart "[ -d /home/bart/PycharmProjects/infinite_cylinder/data/${TOOL}/FerHofSqu1/ ]"
#	then
#		echo
#		echo ">>> Dry run from bart@dart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/"
#		echo
#		rsync -nvtzhre ssh bart@dart:${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/
#	fi
	for PC in ${PCS}
	do
		# shellcheck disable=SC2029
		if ssh ubuntu@"${PC}" "[ -d /home/ubuntu/PycharmProjects/infinite_cylinder/data/${TOOL}/FerHofSqu1/ ]"
		then
			echo
			echo ">>> Dry run from ubuntu@${PC}:${UBUNTU_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/"
			echo
			rsync -nvtzhre ssh ubuntu@"${PC}":${UBUNTU_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/
		fi
	done
done

echo
read -r -p ">>> Are you sure that you want to continue with rsync? [y/N] " response

# actual run
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
    for TOOL in ${TOOLS}
	do
		# shellcheck disable=SC2029
		if ssh bandrews@dirac "[ -d /home/fkp/bandrews/PycharmProjects/infinite_cylinder/data/${TOOL}/FerHofSqu1/ ]"
		then
			echo
			echo ">>> Actual run from bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/"
			echo
			rsync -vtzhre ssh bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/
		fi
		# shellcheck disable=SC2029
#		if ssh bart@dart "[ -d /home/bart/PycharmProjects/infinite_cylinder/data/${TOOL}/FerHofSqu1/ ]"
#		then
#			echo
#			echo ">>> Actual run from bart@dart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/"
#			echo
#			rsync -vtzhre ssh bart@dart:${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/
#		fi
		for PC in ${PCS}
		do
			# shellcheck disable=SC2029
			if ssh ubuntu@"${PC}" "[ -d /home/ubuntu/PycharmProjects/infinite_cylinder/data/${TOOL}/FerHofSqu1/ ]"
			then
				echo
				echo ">>> Actual run from ubuntu@${PC}:${UBUNTU_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/FerHofSqu1/"
				echo
				rsync -vtzhre ssh ubuntu@"${PC}":${UBUNTU_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/FerHofSqu1/
			fi
		done
	done
else
    exit
fi