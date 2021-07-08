# -*- coding: utf-8 -*-
"""
Created on Thu May  7 22:24:36 2020

@author: embog
"""

import picos

# Create a positive semidefinite constant P.
_P = picos.Constant([
    [ 1  -1j,  2  +2j,  1     ],
    [     3j,     -2j, -1  -1j],
    [ 1  +2j, -0.5+1j,  1.5   ]])
P = (_P*_P.H).renamed("P")

# Create a positive semidefinite constant Q.
_Q = picos.Constant([
    [-1  -2j,      2j,  1.5   ],
    [ 1  +2j,     -2j,  2.0-3j],
    [ 1  +2j, -1  +1j,  1  +4j]])
Q = (_Q*_Q.H).renamed("Q")

# Define the problem.
F = picos.Problem()
Z = picos.ComplexVariable("Z", P.shape)
F.set_objective("max", 0.5*picos.trace(Z + Z.H))
F.add_constraint(((P & Z) // (Z.H & Q)) >> 0)

print(F)

# Solve the problem.
F.solve(solver = "mosek", verbosity=2, primals=None)


print("\nOptimal value:", round(F, 4))
print("Optimal Z:", Z.value, sep="\n")
