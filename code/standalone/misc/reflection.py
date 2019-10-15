# A script to reflect the entanglement flow so that only half needs to be calculated

import numpy as np
import csv
import matplotlib.pyplot as plt
from operator import itemgetter
import sys

if __name__ == "__main__":

    input = "haldane_FCI_ent_spec_flow_half.dat"
    output = "haldane_FCI_ent_spec_flow_reflected.dat"

    ent_spec_flow_reflected_file = output
    open(ent_spec_flow_reflected_file, "w")
    ent_spec_flow_reflected_data = open(ent_spec_flow_reflected_file, "a", buffering=1)

    charge = []
    phi = []
    energy = []

    with open(input, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            charge.append(int(row[0]))
            phi.append(float(row[1]))
            energy.append(float(row[2]))
            
    # charge pairs
    pair1 = [0, 1]
    pair2 = [-1, 2]
    pair3 = [-2, 3]
    pair4 = [-3, 4]
    pair5 = [-4, 5]
    pair6 = [-5, 6]
    pair7 = [-6, 7]

    phi_max = 2
    tot_numb_points = 11
    
    for i in range(len(phi)):
        # swap charges
        if charge[i] == pair1[0]:
            charge.append(pair1[1])
        elif charge[i] == pair1[1]:
            charge.append(pair1[0])
        elif charge[i] == pair2[0]:
            charge.append(pair2[1])
        elif charge[i] == pair2[1]:
            charge.append(pair2[0])
        elif charge[i] == pair3[0]:
            charge.append(pair3[1])
        elif charge[i] == pair3[1]:
            charge.append(pair3[0])
        elif charge[i] == pair4[0]:
            charge.append(pair4[1])
        elif charge[i] == pair4[1]:
            charge.append(pair4[0])
        elif charge[i] == pair5[0]:
            charge.append(pair5[1])
        elif charge[i] == pair5[1]:
            charge.append(pair5[0])
        elif charge[i] == pair6[0]:
            charge.append(pair6[1])
        elif charge[i] == pair6[1]:
            charge.append(pair6[0])
        elif charge[i] == pair7[0]:
            charge.append(pair7[1])
        elif charge[i] == pair7[1]:
            charge.append(pair7[0])
        else:
            charge.append(charge[i])

        # reflect phis
        phi.append(float(phi_max-phi[i]))
        energy.append(energy[i])

    ###########################
    # calculate middle column #
    ###########################
    
    energies_left = []
    energies_right = []

    for i in range(len(phi)):
        if int(round(phi[i]*tot_numb_points)) == int(round((tot_numb_points-1)/2-1)):
            energies_left.append(np.array([charge[i], energy[i]]))
        if int(round(phi[i] * tot_numb_points)) == int(round((tot_numb_points-1)/2+1)):
            energies_right.append(np.array([charge[i], energy[i]]))

    sorted_energies_left = sorted(energies_left, key=itemgetter(1))
    sorted_energies_right = sorted(energies_right, key=itemgetter(1))

    for i in range(len(sorted_energies_left)):
        if sorted_energies_left[i][0] == 0:
            for j in range(i+1, len(sorted_energies_left)):
                if sorted_energies_left[j][0] == 1:
                    charge.append(int(sorted_energies_left[i][0]))
                    phi.append(phi_max/2)
                    energy.append((sorted_energies_left[i][1]+sorted_energies_left[j][1])/2)
                    charge.append(int(sorted_energies_left[j][0]))
                    phi.append(phi_max/2)
                    energy.append((sorted_energies_left[i][1] + sorted_energies_left[j][1]) / 2)
                    break
        if int(sorted_energies_left[i][0]) == -1:
            for j in range(i+1, len(sorted_energies_left)):
                if int(sorted_energies_left[j][0]) == 2:
                    charge.append(int(sorted_energies_left[i][0]))
                    phi.append(phi_max/2)
                    energy.append((sorted_energies_left[i][1]+sorted_energies_left[j][1])/2)
                    charge.append(int(sorted_energies_left[j][0]))
                    phi.append(phi_max/2)
                    energy.append((sorted_energies_left[i][1] + sorted_energies_left[j][1]) / 2)
                    break
        if sorted_energies_left[i][0] == -2:
            for j in range(i+1, len(sorted_energies_left)):
                if sorted_energies_left[j][0] == 3:
                    charge.append(int(sorted_energies_left[i][0]))
                    phi.append(phi_max/2)
                    energy.append((sorted_energies_left[i][1]+sorted_energies_left[j][1])/2)
                    charge.append(int(sorted_energies_left[j][0]))
                    phi.append(phi_max/2)
                    energy.append((sorted_energies_left[i][1] + sorted_energies_left[j][1]) / 2)
                    break
        if int(sorted_energies_left[i][0]) == -3:
            for j in range(i+1, len(sorted_energies_left)):
                if int(sorted_energies_left[j][0]) == 4:
                    charge.append(int(sorted_energies_left[i][0]))
                    phi.append(phi_max/2)
                    energy.append((sorted_energies_left[i][1]+sorted_energies_left[j][1])/2)
                    charge.append(int(sorted_energies_left[j][0]))
                    phi.append(phi_max/2)
                    energy.append((sorted_energies_left[i][1] + sorted_energies_left[j][1]) / 2)
                    break

    for i in range(len(phi)):
        # print("{charge:d}\t{phi:.15f}\t{energy:.15f}".format(charge=charge[i], phi=phi[i], energy=energy[i]))
        ent_spec_flow_reflected_data.write("%d\t%.15f\t%.15f\n" % (charge[i], phi[i], energy[i]))