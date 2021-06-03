# --- python imports
import numpy as np


###########################################################################
# interaction_strength (return the interaction potential at a distance r) #
###########################################################################


def interaction_strength(lattice, Vvalue, Vtype_value, rindex):

    if lattice == "Squ":
        r = [1, np.sqrt(2), 2, np.sqrt(5), np.sqrt(8), 3, np.sqrt(10), np.sqrt(13), 4, np.sqrt(17),
             np.sqrt(18), np.sqrt(20), 5, np.sqrt(26), np.sqrt(29), np.sqrt(32), np.sqrt(34), 6, np.sqrt(37), np.sqrt(40),
             np.sqrt(41), np.sqrt(45), 7, np.sqrt(50), np.sqrt(52), np.sqrt(53), np.sqrt(58), np.sqrt(61), 8, np.sqrt(65),
             np.sqrt(68), np.sqrt(72), np.sqrt(73), np.sqrt(74), np.sqrt(80), 9, np.sqrt(82), np.sqrt(85), np.sqrt(89), np.sqrt(90),
             np.sqrt(97), np.sqrt(98), 10, np.sqrt(101), np.sqrt(104), np.sqrt(106), np.sqrt(109), np.sqrt(113), np.sqrt(116), np.sqrt(117)]
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
            neighbors = [(0, 0, np.array([4, 1])), (0, 0, np.array([4, -1])),
                         (0, 0, np.array([1, 4])), (0, 0, np.array([1, -4]))]
        elif n == 11:
            neighbors = [(0, 0, np.array([3, 3])), (0, 0, np.array([3, -3]))]
        elif n == 12:
            neighbors = [(0, 0, np.array([4, 2])), (0, 0, np.array([4, -2])),
                         (0, 0, np.array([2, 4])), (0, 0, np.array([2, -4]))]
        elif n == 13:
            neighbors = [(0, 0, np.array([4, 3])), (0, 0, np.array([4, -3])),
                         (0, 0, np.array([3, 4])), (0, 0, np.array([3, -4])),
                         (0, 0, np.array([5, 0])), (0, 0, np.array([0, 5]))]
        elif n == 14:
            neighbors = [(0, 0, np.array([5, 1])), (0, 0, np.array([5, -1])),
                         (0, 0, np.array([1, 5])), (0, 0, np.array([1, -5]))]
        elif n == 15:
            neighbors = [(0, 0, np.array([5, 2])), (0, 0, np.array([5, -2])),
                         (0, 0, np.array([2, 5])), (0, 0, np.array([2, -5]))]
        elif n == 16:
            neighbors = [(0, 0, np.array([4, 4])), (0, 0, np.array([4, -4]))]
        elif n == 17:
            neighbors = [(0, 0, np.array([5, 3])), (0, 0, np.array([5, -3])),
                         (0, 0, np.array([3, 5])), (0, 0, np.array([3, -5]))]
        elif n == 18:
            neighbors = [(0, 0, np.array([6, 0])), (0, 0, np.array([0, 6]))]
        elif n == 19:
            neighbors = [(0, 0, np.array([6, 1])), (0, 0, np.array([6, -1])),
                         (0, 0, np.array([1, 6])), (0, 0, np.array([1, -6]))]
        elif n == 20:
            neighbors = [(0, 0, np.array([6, 2])), (0, 0, np.array([6, -2])),
                         (0, 0, np.array([2, 6])), (0, 0, np.array([2, -6]))]
        elif n == 21:
            neighbors = [(0, 0, np.array([5, 4])), (0, 0, np.array([5, -4])),
                         (0, 0, np.array([4, 5])), (0, 0, np.array([4, -5]))]
        elif n == 22:
            neighbors = [(0, 0, np.array([6, 3])), (0, 0, np.array([6, -3])),
                         (0, 0, np.array([3, 6])), (0, 0, np.array([3, -6]))]
        elif n == 23:
            neighbors = [(0, 0, np.array([7, 0])), (0, 0, np.array([0, 7]))]
        elif n == 24:
            neighbors = [(0, 0, np.array([5, 5])), (0, 0, np.array([5, -5])),
                         (0, 0, np.array([7, 1])), (0, 0, np.array([7, -1])),
                         (0, 0, np.array([1, 7])), (0, 0, np.array([1, -7]))]
        elif n == 25:
            neighbors = [(0, 0, np.array([6, 4])), (0, 0, np.array([6, -4])),
                         (0, 0, np.array([4, 6])), (0, 0, np.array([4, -6]))]
        elif n == 26:
            neighbors = [(0, 0, np.array([7, 2])), (0, 0, np.array([7, -2])),
                         (0, 0, np.array([2, 7])), (0, 0, np.array([2, -7]))]
        elif n == 27:
            neighbors = [(0, 0, np.array([7, 3])), (0, 0, np.array([7, -3])),
                         (0, 0, np.array([3, 7])), (0, 0, np.array([3, -7]))]
        elif n == 28:
            neighbors = [(0, 0, np.array([6, 5])), (0, 0, np.array([6, -5])),
                         (0, 0, np.array([5, 6])), (0, 0, np.array([5, -6]))]
        elif n == 29:
            neighbors = [(0, 0, np.array([8, 0])), (0, 0, np.array([0, 8])),
                         (0, 0, np.array([7, 4])), (0, 0, np.array([7, -4])),
                         (0, 0, np.array([4, 7])), (0, 0, np.array([4, -7]))]
        elif n == 30:
            neighbors = [(0, 0, np.array([8, 1])), (0, 0, np.array([8, -1])),
                         (0, 0, np.array([1, 8])), (0, 0, np.array([1, -8]))]
        elif n == 31:
            neighbors = [(0, 0, np.array([8, 2])), (0, 0, np.array([8, -2])),
                         (0, 0, np.array([2, 8])), (0, 0, np.array([2, -8]))]
        elif n == 32:
            neighbors = [(0, 0, np.array([6, 6])), (0, 0, np.array([6, -6]))]
        elif n == 33:
            neighbors = [(0, 0, np.array([8, 3])), (0, 0, np.array([8, -3])),
                         (0, 0, np.array([3, 8])), (0, 0, np.array([3, -8]))]
        elif n == 34:
            neighbors = [(0, 0, np.array([7, 5])), (0, 0, np.array([7, -5])),
                         (0, 0, np.array([5, 7])), (0, 0, np.array([5, -7]))]
        elif n == 35:
            neighbors = [(0, 0, np.array([8, 4])), (0, 0, np.array([8, -4])),
                         (0, 0, np.array([4, 8])), (0, 0, np.array([4, -8]))]
        elif n == 36:
            neighbors = [(0, 0, np.array([9, 0])), (0, 0, np.array([0, 9]))]
        elif n == 37:
            neighbors = [(0, 0, np.array([9, 1])), (0, 0, np.array([9, -1])),
                         (0, 0, np.array([1, 9])), (0, 0, np.array([1, -9]))]
        elif n == 38:
            neighbors = [(0, 0, np.array([7, 6])), (0, 0, np.array([7, -6])),
                         (0, 0, np.array([6, 7])), (0, 0, np.array([6, -7])),
                         (0, 0, np.array([9, 2])), (0, 0, np.array([9, -2])),
                         (0, 0, np.array([2, 9])), (0, 0, np.array([2, -9]))]
        elif n == 39:
            neighbors = [(0, 0, np.array([8, 5])), (0, 0, np.array([8, -5])),
                         (0, 0, np.array([5, 8])), (0, 0, np.array([5, -8]))]
        elif n == 40:
            neighbors = [(0, 0, np.array([9, 3])), (0, 0, np.array([9, -3])),
                         (0, 0, np.array([3, 9])), (0, 0, np.array([3, -9]))]
        elif n == 41:
            neighbors = [(0, 0, np.array([9, 4])), (0, 0, np.array([9, -4])),
                         (0, 0, np.array([4, 9])), (0, 0, np.array([4, -9]))]
        elif n == 42:
            neighbors = [(0, 0, np.array([7, 7])), (0, 0, np.array([7, -7]))]
        elif n == 43:
            neighbors = [(0, 0, np.array([8, 6])), (0, 0, np.array([8, -6])),
                         (0, 0, np.array([6, 8])), (0, 0, np.array([6, -8])),
                         (0, 0, np.array([10, 0])), (0, 0, np.array([0, 10]))]
        elif n == 44:
            neighbors = [(0, 0, np.array([10, 1])), (0, 0, np.array([10, -1])),
                         (0, 0, np.array([1, 10])), (0, 0, np.array([1, -10]))]
        elif n == 45:
            neighbors = [(0, 0, np.array([10, 2])), (0, 0, np.array([10, -2])),
                         (0, 0, np.array([2, 10])), (0, 0, np.array([2, -10]))]
        elif n == 46:
            neighbors = [(0, 0, np.array([9, 5])), (0, 0, np.array([9, -5])),
                         (0, 0, np.array([5, 9])), (0, 0, np.array([5, -9]))]
        elif n == 47:
            neighbors = [(0, 0, np.array([10, 3])), (0, 0, np.array([10, -3])),
                         (0, 0, np.array([3, 10])), (0, 0, np.array([3, -10]))]
        elif n == 48:
            neighbors = [(0, 0, np.array([8, 7])), (0, 0, np.array([8, -7])),
                         (0, 0, np.array([7, 8])), (0, 0, np.array([7, -8]))]
        elif n == 49:
            neighbors = [(0, 0, np.array([10, 4])), (0, 0, np.array([10, -4])),
                         (0, 0, np.array([4, 10])), (0, 0, np.array([4, -10]))]
        elif n == 50:
            neighbors = [(0, 0, np.array([9, 6])), (0, 0, np.array([9, -6])),
                         (0, 0, np.array([6, 9])), (0, 0, np.array([6, -9]))]
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
