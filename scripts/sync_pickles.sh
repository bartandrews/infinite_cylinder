#!/bin/bash

# This script will sync all pickles to dart

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
PCS="baandr1 baandr2 baandr3"  # username == ubuntu
SRC_BASE=/home/ubuntu
DST_BASE=/home/bart
DIR_PATH=PycharmProjects/infinite_cylinder/pickles/ground_state

# dry run
for DIR in ${MODELS}
do
	for PC in ${PCS}
	do
		if ssh ubuntu@${PC} "[ -d /home/ubuntu/PycharmProjects/infinite_cylinder/pickles/ground_state/${DIR}/ ]"
		then
			echo
			echo ">>> Dry run from ubuntu@${PC}:${SRC_BASE}/${DIR_PATH}/${DIR}/ to bart@dart:${DST_BASE}/${DIR_PATH}/${DIR}/"
			echo
			ssh -A dart rsync -nvtzhre ssh ubuntu@${PC}:${SRC_BASE}/${DIR_PATH}/${DIR}/ ${DST_BASE}/${DIR_PATH}/${DIR}/
		fi
	done
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
			if ssh ubuntu@${PC} "[ -d /home/ubuntu/PycharmProjects/infinite_cylinder/pickles/ground_state/${DIR}/ ]"
			then
				echo
				echo ">>> Actual run from ubuntu@${PC}:${SRC_BASE}/${DIR_PATH}/${DIR}/ to bart@dart:${DST_BASE}/${DIR_PATH}/${DIR}/"
				echo
				ssh -A dart rsync -vtzhre ssh ubuntu@${PC}:${SRC_BASE}/${DIR_PATH}/${DIR}/ ${DST_BASE}/${DIR_PATH}/${DIR}/
			fi
		done
	done
else
    exit
fi