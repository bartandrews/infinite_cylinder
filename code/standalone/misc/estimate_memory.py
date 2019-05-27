#!/usr/bin/env python3

import numpy as np

#  D = [22, 28, 22, 28, 22, 28, 22, 28, 22]
N = 12  # number of sites in the unit cell
D = [68]*N   # what M.H_MPO.chi gives
chi = 400 # max_chi in trunc_params
d = 16  # local dimension of the Site

psi = chi**2 * d * (N + 20)
env = chi**2 * (sum(D) + 2 * max(D))
print("total # floats", psi + env)
print("GB for complex:", (psi + env) / (1024**3) * 128/8)
