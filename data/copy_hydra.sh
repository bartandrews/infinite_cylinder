#!/bin/bash

for tool in charge_pump corr_len ent_scal ent_spec_flow ent_spec_mom ent_spec_real ent_spec_V_flow
do
    scp s3it:~/infinite_cylinder/data/$tool/*.dat $tool
done

