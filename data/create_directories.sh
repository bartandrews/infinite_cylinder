#!/bin/bash

for tool in charge_pump corr_len_kappa_flow double_occ ent_spec_flow ent_spec_mom ent_spec_V_flow corr_len corr_len_tau_flow energy ent_spec_invt2dash_flow ent_spec_real overlap corr_len_invt2dash_flow ent_scal ent_spec_kappa_flow ent_spec_tau_flow
do
    for model in BosonicHofstadter FermionicHofstadter BosonicHex1 FermionicHex1 BosonicHex1Hex5 FermionicHex1Hex5 BosonicHex1Hex5Orbital FermionicHex1Hex5Orbital 
    do
        #mkdir -p $tool/$model
        printf '%s\n%s' '*' '!.gitignore' > $tool/$model/.gitignore
    done
done
