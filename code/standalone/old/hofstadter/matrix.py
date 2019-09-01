import numpy
from numpy import exp, zeros, pi, cos


class PerturbedMatrix(object):
    def __init__(self, g, a, size):
        super(PerturbedMatrix, self).__init__()
        self.matrix = construct_matrix(g, a, size)
        self.alphas = [a for i in range(0, size)]

    @property
    def eigenvalues(self):
        return numpy.linalg.eigvals(self.matrix)


def construct_matrix(g, a, size):
    _matrix = zeros(shape=(size, size))
    matrix_norm = (2 * pi * a)
    minus_computed = exp(-g)
    plus_computed = exp(g)

    for i in range(0, size):
        _matrix[i, i] = 2*cos(matrix_norm*i)
        minus, plus = _construct_indeces(i, size)

        if minus is not None:
            _matrix[i, minus] = minus_computed

        if plus is not None:
            _matrix[i, plus] = plus_computed

    return _matrix

def _calculate_eigenvalues(matrix):
    return numpy.linalg.eigvals(matrix)

def matrix_eigenvalues(g, a, size):
    matrix = construct_matrix(g, a, size)

    print(matrix)

    return _calculate_eigenvalues(matrix)


def _construct_indeces(index, size):
    if index == 0:
        return size - 1, 1
    if index == size - 1:
        return index - 1, 0

    return index - 1, index + 1