# A script to rotate the ent_spec_mom so that the edge states are visible

import numpy as np
import csv
import matplotlib.pyplot as plt
from operator import itemgetter
import sys

if __name__ == "__main__":

    input = "/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_mom/FerHofSqu1/ent_spec_mom_FerHofSqu1_chi_500_chiK_500_t1_1_V_3_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat"
    output = "/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_mom/FerHofSqu1/ent_spec_mom_FerHofSqu1_chi_500_chiK_500_t1_1_V_3_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat.rotated"

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
    Ly = 10
    # desired rotation to the left (integer)
    shift = 6

    for i in range(len(mom)):

        # the units are already pi
        mom[i] -= shift*2/Ly
        if mom[i] < -1:
            mom[i] += 2

        ent_spec_mom_rotated_data.write("%d\t%.15f\t%.15f\n" % (charge[i], mom[i], energy[i]))
