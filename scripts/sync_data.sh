#!/bin/bash
# sync_data.sh --- sync data from various to laptop (execute on laptop)

# rsync options
#
# n = dry run (optional)
# c = checksum
# v = verbose
# t = preserve modification time
# z = compress files for transfer
# h = human-readable output
# r = recursive (optional)
# e ssh = transfer files over a secure connection
# ignore-existing = skip updating files that exist on receiver

MODELS="BosHofSqu1 FerHofSqu1"  # the model directories that we would like to sync
TOOLS="ent_spec_real ent_spec_mom density corr_func corr_func_ext overlap charge_pump ent_spec_flow ent_spec_V_flow ent_V_flow corr_len_V_flow ent_corr_len energy_V_flow"  # defined in code
PCS="baandr1 baandr2 baandr3"  # remote hostnames that we intend to run on (username == ubuntu is assumed throughout)
BANDREWS_BASE=/home/fkp/bandrews  # the home directory on the remote computers with username == bandrews
BART_BASE=/home/bart  # the home directory on the local computer with username == bart
UBUNTU_BASE=/home/ubuntu  # the home directory on the remote computers with username == ubuntu
DIR_PATH_HOME=PycharmProjects/infinite_cylinder/data  # path to the ground_state directory relative to home

# dry run
for DIR in ${MODELS}
do
  for TOOL in ${TOOLS}
  do
#    for PROJECT in project1 project2
#    do
#      # shellcheck disable=SC2029
#      if ssh bandrews@dirac "[ -d /disk/data11/tfp/BartMadhav/${PROJECT}/data/${TOOL}/${DIR}/ ]"
#      then
#        echo
#        echo ">>> Dry run from bandrews@dirac:/disk/data11/tfp/BartMadhav/${PROJECT}/data/${TOOL}/${DIR}/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/"
#        echo
#        rsync -cnvtzhre ssh --ignore-existing bandrews@dirac:/disk/data11/tfp/BartMadhav/${PROJECT}/data/"${TOOL}"/"${DIR}"/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/
#      fi
#    done
    # shellcheck disable=SC2029
#    if ssh bandrews@dirac "[ -d /home/fkp/bandrews/PycharmProjects/infinite_cylinder/data/${TOOL}/${DIR}/ ]"
#    then
#      echo
#      echo ">>> Dry run from bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/"
#      echo
#      rsync -cnvtzhre ssh --ignore-existing bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/
#    fi
    # shellcheck disable=SC2029
#  	if ssh bart@dart "[ -d /home/bart/PycharmProjects/infinite_cylinder/data/${TOOL}/${DIR}/ ]"
#  	then
#  		echo
#  		echo ">>> Dry run from bart@dart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/"
#  		echo
#  		rsync -cnvtzhre ssh --ignore-existing bart@dart:${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/
#  	fi
    for PC in ${PCS}
    do
      # shellcheck disable=SC2029
      if ssh ubuntu@"${PC}" "[ -d /home/ubuntu/PycharmProjects/infinite_cylinder/data/${TOOL}/${DIR}/ ]"
      then
        echo
        echo ">>> Dry run from ubuntu@${PC}:${UBUNTU_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/"
        echo
        rsync -cnvtzhre ssh --ignore-existing ubuntu@"${PC}":${UBUNTU_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/
      fi
    done
  done
done

echo
read -r -p ">>> Are you sure that you want to continue with rsync? [y/N] " response

# actual run
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
  for DIR in ${MODELS}
  do
    for TOOL in ${TOOLS}
    do
#      for PROJECT in project1 project2
#      do
#        # shellcheck disable=SC2029
#        if ssh bandrews@dirac "[ -d /disk/data11/tfp/BartMadhav/${PROJECT}/data/${TOOL}/${DIR}/ ]"
#        then
#          echo
#          echo ">>> Dry run from bandrews@dirac:/disk/data11/tfp/BartMadhav/${PROJECT}/data/${TOOL}/${DIR}/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/"
#          echo
#          rsync -cvtzhre ssh --ignore-existing bandrews@dirac:/disk/data11/tfp/BartMadhav/${PROJECT}/data/"${TOOL}"/"${DIR}"/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/
#        fi
#      done
      # shellcheck disable=SC2029
#      if ssh bandrews@dirac "[ -d /home/fkp/bandrews/PycharmProjects/infinite_cylinder/data/${TOOL}/${DIR}/ ]"
#      then
#        echo
#        echo ">>> Actual run from bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/"
#        echo
#        rsync -cvtzhre ssh --ignore-existing bandrews@dirac:${BANDREWS_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/
#      fi
      # shellcheck disable=SC2029
#      if ssh bart@dart "[ -d /home/bart/PycharmProjects/infinite_cylinder/data/${TOOL}/${DIR}/ ]"
#      then
#        echo
#        echo ">>> Actual run from bart@dart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/"
#        echo
#        rsync -cvtzhre ssh --ignore-existing bart@dart:${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/
#      fi
      for PC in ${PCS}
      do
        # shellcheck disable=SC2029
        if ssh ubuntu@"${PC}" "[ -d /home/ubuntu/PycharmProjects/infinite_cylinder/data/${TOOL}/${DIR}/ ]"
        then
          echo
          echo ">>> Actual run from ubuntu@${PC}:${UBUNTU_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/ to bart@bart:${BART_BASE}/${DIR_PATH_HOME}/${TOOL}/${DIR}/"
          echo
          rsync -cvtzhre ssh --ignore-existing ubuntu@"${PC}":${UBUNTU_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/ ${BART_BASE}/${DIR_PATH_HOME}/"${TOOL}"/"${DIR}"/
        fi
      done
    done
	done
else
    exit
fi