# --- python imports
import numpy as np


###########################################################################
# interaction_strength (return the interaction potential at a distance r) #
###########################################################################


def interaction_strength(lattice, Vvalue, Vtype_value, rindex):

    if lattice == "Squ":
        r = [1, np.sqrt(2), 2, np.sqrt(5), np.sqrt(8), 3, np.sqrt(10), np.sqrt(13), 4, np.sqrt(18)]
    elif lattice == "Hex":
        r = [1, np.sqrt(3), 2, np.sqrt(7), 3, np.sqrt(12), np.sqrt(13), 4, np.sqrt(19), np.sqrt(21)]
    elif lattice == "Tri":
        r = [1, np.sqrt(3), 2]
    else:
        raise ValueError("Unknown lattice parameter.")

    if Vtype_value == 'Coulomb':
        potential = Vvalue / r[rindex]
    elif Vtype_value == 'Yukawa':
        potential = np.exp(-r[rindex]) * Vvalue / r[rindex]
    else:
        raise ValueError("Unknown Vtype for the interaction_strength function.")
    return potential


##################################################################
# NN (list of neighbors for a given lattice at a specific range) #
##################################################################


def NN(lattice, n):
    if lattice == "Squ":
        if n == 1:
            neighbors = [(0, 0, np.array([1, 0])), (0, 0, np.array([0, 1]))]
        elif n == 2:
            neighbors = [(0, 0, np.array([1, 1])), (0, 0, np.array([1, -1]))]
        elif n == 3:
            neighbors = [(0, 0, np.array([2, 0])), (0, 0, np.array([0, 2]))]
        elif n == 4:
            neighbors = [(0, 0, np.array([2, 1])), (0, 0, np.array([2, -1])),
                         (0, 0, np.array([1, 2])), (0, 0, np.array([1, -2]))]
        elif n == 5:
            neighbors = [(0, 0, np.array([2, 2])), (0, 0, np.array([2, -2]))]
        elif n == 6:
            neighbors = [(0, 0, np.array([3, 0])), (0, 0, np.array([0, 3]))]
        elif n == 7:
            neighbors = [(0, 0, np.array([3, 1])), (0, 0, np.array([3, -1])),
                         (0, 0, np.array([1, 3])), (0, 0, np.array([1, -3]))]
        elif n == 8:
            neighbors = [(0, 0, np.array([3, 2])), (0, 0, np.array([3, -2])),
                         (0, 0, np.array([2, 3])), (0, 0, np.array([2, -3]))]
        elif n == 9:
            neighbors = [(0, 0, np.array([4, 0])), (0, 0, np.array([0, 4]))]
        elif n == 10:
            neighbors = [(0, 0, np.array([3, 3])), (0, 0, np.array([3, -3]))]
        else:
            raise ValueError("nearest_neighbors not implemented for this value of n.")
    elif lattice == "Hex":
        if n == 1:
            neighbors = [(0, 1, np.array([0, 0])), (1, 0, np.array([1, 0])), (1, 0, np.array([0, 1]))]
        elif n == 2:
            neighbors = [(0, 0, np.array([1, 0])), (0, 0, np.array([0, 1])), (0, 0, np.array([1, -1])),
                         (1, 1, np.array([1, 0])), (1, 1, np.array([0, 1])), (1, 1, np.array([1, -1]))]
        elif n == 3:
            neighbors = [(1, 0, np.array([1, 1])), (0, 1, np.array([-1, 1])), (0, 1, np.array([1, -1]))]
        elif n == 4:
            neighbors = [(0, 1, np.array([0, 1])), (0, 1, np.array([1, 0])),
                         (0, 1, np.array([1, -2])), (0, 1, np.array([0, -2])),
                         (0, 1, np.array([-2, 0])), (0, 1, np.array([-2, 1]))]
        elif n == 5:
            neighbors = [(0, 0, np.array([1, 1])), (0, 0, np.array([2, -1])), (0, 0, np.array([-1, 2])),
                         (1, 1, np.array([1, 1])), (1, 1, np.array([2, -1])), (1, 1, np.array([-1, 2]))]
        elif n == 6:
            neighbors = [(0, 0, np.array([2, 0])), (0, 0, np.array([0, 2])), (0, 0, np.array([2, -2])),
                         (1, 1, np.array([2, 0])), (1, 1, np.array([0, 2])), (1, 1, np.array([2, -2]))]
        elif n == 7:
            neighbors = [(0, 1, np.array([-1, 2])), (0, 1, np.array([-2, 2])),
                         (0, 1, np.array([2, -1])), (0, 1, np.array([2, -2])),
                         (0, 1, np.array([-2, -1])), (0, 1, np.array([-1, -2]))]
        elif n == 8:
            neighbors = [(0, 1, np.array([1, 1])), (0, 1, np.array([1, -3])), (0, 1, np.array([-3, 1]))]
        elif n == 9:
            neighbors = [(0, 1, np.array([0, 2])), (0, 1, np.array([2, 0])),
                         (0, 1, np.array([2, -3])), (0, 1, np.array([0, -3])),
                         (0, 1, np.array([-3, 0])), (0, 1, np.array([-3, 2]))]
        elif n == 10:
            neighbors = [(0, 0, np.array([1, 2])), (0, 0, np.array([2, 1])), (0, 0, np.array([3, -1])),
                         (0, 0, np.array([3, -2])), (0, 0, np.array([2, -3])), (0, 0, np.array([1, -3])),
                         (1, 1, np.array([1, 2])), (1, 1, np.array([2, 1])), (1, 1, np.array([3, -1])),
                         (1, 1, np.array([3, -2])), (1, 1, np.array([2, -3])), (1, 1, np.array([1, -3]))]
        else:
            raise ValueError("nearest_neighbors not implemented for this value of n.")
    elif lattice == "Tri":
        if n == 1:
            neighbors = [(0, 0, np.array([1, 0])), (0, 0, np.array([-1, 1])), (0, 0, np.array([0, -1]))]
        elif n == 2:
            neighbors = [(0, 0, np.array([1, 1])), (0, 0, np.array([-2, 1])), (0, 0, np.array([1, -2]))]
        elif n == 3:
            neighbors = [(0, 0, np.array([2, 0])), (0, 0, np.array([-2, 2])), (0, 0, np.array([0, -2]))]
        else:
            raise ValueError("nearest_neighbors not implemented for this value of n.")
    else:
        raise ValueError("Unable to infer lattice from model name.")

    return neighbors
