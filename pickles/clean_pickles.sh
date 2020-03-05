#!/bin/bash

for flow in phi_flow Ly_flow U_flow V_flow kappa_flow
do
    find $flow/ -type f -size 0 -delete
done

