#!/bin/bash

read -r -p "Are you sure that you want to delete all files in data, logs, and pickles? [y/N] " response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
  cd data
  for tool in overlap charge_pump ent_spec_flow ent_scal ent_spec_real ent_spec_mom corr_len_U_flow double_occ_U_flow corr_len_V_flow ent_spec_V_flow corr_len_kappa_flow ent_spec_kappa_flow
  do
      rm -rf $tool/
  done

  cd logs
  for flow in phi_flow Ly_flow U_flow V_flow kappa_flow
  do
      rm -rf $flow/
  done

  cd pickles
  for flow in phi_flow Ly_flow U_flow V_flow kappa_flow
  do
      rm -rf $flow/
  done
else
    exit
fi







