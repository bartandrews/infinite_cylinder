#!/bin/bash

#for tool in charge_pump corr_len double_occ energy ent_scal ent_spec_flow ent_spec_mom ent_spec_real ent_spec_V_flow overlap 
for tool in corr_len_invt2dash_flow ent_spec_invt2dash_flow 
do
    for model in BosonicHofstadter FermionicHofstadter BosonicHex1 FermionicHex1 BosonicHex1Hex5 FermionicHex1Hex5 BosonicHex1Hex5Orbital FermionicHex1Hex5Orbital 
    do
        mkdir -p $tool/$model
        printf '%s\n%s' '*' '!.gitignore' > $tool/$model/.gitignore
    done
done
