# A script to rotate the ent_spec_mom so that the edge states are visible

import numpy as np
import csv
import matplotlib.pyplot as plt
from operator import itemgetter
import sys

if __name__ == "__main__":

    input = "../main_figures_rev1/Ly_flow_mag_detail/ent_spec_mom_Ly_9.dat"
    output = "../main_figures_rev1/Ly_flow_mag_detail/ent_spec_mom_Ly_9.dat.rotated"

    ent_spec_mom_rotated_file = output
    open(ent_spec_mom_rotated_file, "w")
    ent_spec_mom_rotated_data = open(ent_spec_mom_rotated_file, "a", buffering=1)

    charge = []
    mom = []
    energy = []

    with open(input, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            charge.append(int(row[0]))
            mom.append(float(row[1]))
            energy.append(float(row[2]))

    # original Ly of file
    Ly = 9
    # desired rotation to the left (integer)
    shift = 7

    for i in range(len(mom)):

        # the units are already pi
        mom[i]-=shift*2/Ly
        if mom[i]<-1:
            mom[i]+=2

        ent_spec_mom_rotated_data.write("%d\t%.15f\t%.15f\n" % (charge[i], mom[i], energy[i]))