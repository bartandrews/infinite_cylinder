#!/bin/bash

for tool in charge_pump corr_len double_occ energy ent_scal ent_spec_flow ent_spec_mom ent_spec_real ent_spec_V_flow overlap 
do
    for model in BosonicHofstadter FermionicHofstadter BosonicHex1 FermionicHex1 BosonicHex1Hex5 FermionicHex1Hex5 BosonicHex1Hex5Orbital FermionicHex1Hex5Orbital 
    do
	for tag in nphi_1_4 nphi_8_9 nphi_10_11 
        do
            rm -f $tool/$model/*.$tag
        done
    done
done
