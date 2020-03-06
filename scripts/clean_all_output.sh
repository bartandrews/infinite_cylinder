#!/bin/bash

cd ../data
for tool in overlap charge_pump ent_spec_flow ent_scal ent_spec_real ent_spec_mom corr_len_U_flow double_occ_U_flow corr_len_V_flow ent_spec_V_flow corr_len_kappa_flow ent_spec_kappa_flow
do
    find $tool/ -type f -size 0 -delete
done

cd ../logs
for flow in phi_flow Ly_flow U_flow V_flow kappa_flow
do
    find $flow/ -type f -size 0 -delete
done

cd ../pickles
for flow in phi_flow Ly_flow U_flow V_flow kappa_flow
do
    find $flow/ -type f -size 0 -delete
done

