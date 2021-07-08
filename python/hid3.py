# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:04:28 2020

@author: embog
"""


import picos
import numpy as np
import itertools

OUTPUT_VALUES  = [[0,1],[0,1],[0,1]]
INPUT_VALUES = [[0,1,2],[0,1],[0,1]]
NR_PARTIES = len(OUTPUT_VALUES)

def test_party_nr():
    assert len(INPUT_VALUES) == NR_PARTIES

    
#%%
    
'''
Quantum distribution
'''

import sympy as sym

from sympy.physics.quantum import TensorProduct, Dagger

def norm_pure_state(state):
    # [0] in the output because it outputs a 1x1 matrix annoyingly
    return sym.sqrt((Dagger(state)*state)[0]) 

def projector_onto_state(state):
    normm = norm_pure_state(state)
    if normm != 0:    
        return state*Dagger(state)/(normm**2)  # |psi><psi|/<psi|psi>
    else:
        print("ERROR: norm 0")
        return 0
    
sig0   = sym.Matrix([[1, 0],[0, 1]])
# All these are written in the Z basis, |0>, |1>
sig1   = sym.Matrix([[0, 1],[1, 0]]) 
sig2   = sym.Matrix([[0,-sym.I],[sym.I, 0]])
sig3   = sym.Matrix([[1, 0],[0,-1]])

def ret_proj_eigs(obs):
    obs_vec, obs_val = obs.diagonalize()
    ret = []
    for i in range(obs_vec.shape[1]):
        ret.append(projector_onto_state(obs_vec[:,i]))
    return ret

#Alice
A0 = sig3 #sig_z
A1 = sig1 #sig_x
A2 = sig2 #sig_y

A_obs = [A0,A1,A2]
A_proj = []
for i in range(len(INPUT_VALUES[0])):
    A_proj.append(ret_proj_eigs(A_obs[i]))


#Bob
phi = np.arctan(1/np.sqrt(2))
B0 = np.cos(phi) * sig1 + np.sin(phi) * sig2
B1 = np.cos(phi) * sig1 - np.sin(phi) * sig2

B_obs = [B0,B1]
B_proj = []
for i in range(len(INPUT_VALUES[1])):
    B_proj.append(ret_proj_eigs(B_obs[i]))


#Charie
C0 = sig3
C1 = sig1

C_obs = [C0,C1]
C_proj = []
for i in range(len(INPUT_VALUES[2])):
    C_proj.append(ret_proj_eigs(C_obs[i]))

P_proj = [A_proj,B_proj,C_proj]

#U
e1 = sym.Matrix([1,0])
e2 = sym.Matrix([0,1])



Phi_minus = 1/np.sqrt(2) * (TensorProduct(e1,e1)-TensorProduct(e2,e2))
Psi_plus  = 1/np.sqrt(2) * (TensorProduct(e1,e2)+TensorProduct(e2,e1))
Phi_plus  = 1/np.sqrt(2) * (TensorProduct(e1,e1)+TensorProduct(e2,e2))

ALPHA = sym.symbols("alpha")
ini_state = Phi_plus * Dagger(Phi_plus)
ini_state = (ALPHA) * ini_state + (1-ALPHA) * 1.0/4 * TensorProduct(sig0,sig0)

alpha = np.pi/8
psi0 =  np.sin(alpha) * Phi_minus + np.cos(alpha) * Psi_plus
psi1 = -np.cos(alpha) * Phi_minus + np.sin(alpha) * Psi_plus

# A and B+C share Phi_plus. This U acts on B+C's part. If the box measures 0,
# then it gives psi0 to Bob and Charlie to measure; if the box measures 1, it 
# gives psi1 to B+c.
U = psi0 * Dagger(e1) + psi1 * Dagger(e2)

# Now I will get the final state for all three parties.

Big_U = TensorProduct(sig0,U)
final_state = Big_U * ini_state * Dagger(Big_U)


def meas_proj(party,i,o):
    return P_proj[party][i][o]

def probability(rho, inputs, outputs):
    kronn = meas_proj(0,inputs[0],outputs[0])
    for i in np.arange(1,NR_PARTIES):
        kronn = TensorProduct(kronn,meas_proj(i,inputs[i],outputs[i]))
    ret = sym.simplify(sym.trace(rho*kronn))
    return sym.nsimplify(ret,tolerance=1e-10,rational=False).evalf()

def test_probability():
    for x,y,z in itertools.product(*INPUT_VALUES):
        suma = 0
        for a,b,c in itertools.product(*OUTPUT_VALUES):
            suma = suma + probability(final_state*Dagger(final_state),
                                      [x,y,z],
                                      [a,b,c])
        assert (suma - 1) < 1e-10
        
#%%

nr_det_points = 8
det = [[[0 for a in range(len(OUTPUT_VALUES[0]))]
             for x in range(len(INPUT_VALUES[0]))]
                for i in range(nr_det_points)]

counter = 0
for a_tuple in itertools.product(OUTPUT_VALUES[0],repeat=len(INPUT_VALUES[0])):
    a_list = list(a_tuple)
    for x in INPUT_VALUES[0]:
        det[counter][x][a_tuple[x]] = 1
    counter = counter + 1


P = picos.Problem()


#%%
# CONSTRAINTS
positivity_constraints = []


alpha = picos.RealVariable("alpha")
positivity_constraints = [*positivity_constraints, alpha >= 0, alpha <= 1]


q = picos.RealVariable("q",nr_det_points*2*2*2*2)
positivity_constraints = [*positivity_constraints, q>=0]


P.add_list_of_constraints(positivity_constraints)


q_vars = [[[[[0 for c in range(len(OUTPUT_VALUES[2]))]
                for b in range(len(OUTPUT_VALUES[1]))]
                    for z in range(len(INPUT_VALUES[2]))]
                        for y in range(len(INPUT_VALUES[1]))]
                            for lam in range(nr_det_points)]
i=0
for lam in range(nr_det_points):
    for y in range(len(INPUT_VALUES[1])):
        for z in range(len(INPUT_VALUES[2])):
            for b in range(len(OUTPUT_VALUES[1])):
                for c in range(len(OUTPUT_VALUES[2])):         
                    q_vars[lam][y][z][b][c] = q[i]
                    i = i + 1
    
non_signalling_constraints=[]
        
# non signaling for bob: sum over C, make this not depent on C's input
for lam in range(nr_det_points):
    for b,y in itertools.product(OUTPUT_VALUES[1],INPUT_VALUES[1]):
        for z1,z2 in itertools.combinations(INPUT_VALUES[2],2):
            # combinations('ABCD', 2) --> AB AC AD BC BD CD
            # we use combinations becausee for no signalling we want for
            # every pair of different y1's to show that they are different, 
            # but y1_1 = y1_2 is the same as y1_2 = y1_1, so we use 
            # combinations instead of iter.product to avoid adding the same 
            # equality twice, or redundant equalities such as y1_1 = y1_1.
            
            # we calculate the two marginals for bob by summing over charlie
            suma_z1 = 0
            for c in OUTPUT_VALUES[2]:
                suma_z1 = suma_z1 + q_vars[lam][y][z1][b][c]
            suma_z2 = 0
            for c in OUTPUT_VALUES[2]:
                suma_z2 = suma_z2 + q_vars[lam][y][z2][b][c]  
            non_signalling_constraints.append(suma_z1 == suma_z2)
# non signaling for charlie: sum over B, make this not depent on B's input
for lam in range(nr_det_points):
    for c,z in itertools.product(OUTPUT_VALUES[2],INPUT_VALUES[2]):
        for y1,y2 in itertools.combinations(INPUT_VALUES[1],2):
            suma_y1 = 0
            for b in OUTPUT_VALUES[1]:
                suma_y1 = suma_y1 + q_vars[lam][y1][z][b][c]
            suma_y2 = 0
            for b in OUTPUT_VALUES[1]:
                suma_y2 = suma_y2 + q_vars[lam][y2][z][b][c]
            non_signalling_constraints.append(suma_y1 == suma_y2)
            
P.add_list_of_constraints(non_signalling_constraints)


prob_constraints=[]
prob_coordinates=[]

for a,b,c in itertools.product(*OUTPUT_VALUES):
    for x,y,z in itertools.product(*INPUT_VALUES):
        
        suma = 0
        for i in range(nr_det_points):
            suma = suma + det[i][x][a] * q_vars[i][y][z][b][c]
        
        prob1 = probability(final_state,[x,y,z],[a,b,c])
        prob_coordinates.append([[x,y,z],[a,b,c]])
                
        # change symbolic ALPHA to PICOS variable alpha
        prob2 = (sym.lambdify(ALPHA,prob1))(alpha)

        prob_constraints.append(suma == prob2)

P.add_list_of_constraints(prob_constraints)


#%%
# Now solve the problem.

#P.set_objective("max",alpha)
P.set_objective(None)


solution = P.solve(solver="cvxopt")
print("Problem status:", solution.problemStatus)

print("Noise tolerance: alpha=", alpha.value)

solution
#%%

'''
This is to compare with the things from Joe's notes.
'''

correlator = ( -2*TensorProduct(A2,B0,sig0)+2*TensorProduct(A2,B1,sig0)
              + TensorProduct(A0,B0,C0) + TensorProduct(A0,B1,C1)
              + TensorProduct(A1,B1,C1) - TensorProduct(A1,B0,C0)
              + TensorProduct(A0,B0,C1) + TensorProduct(A0,B1,C0)
              + TensorProduct(A1,B0,C1) - TensorProduct(A1,B1,C0) )

def corr_noise_resistance(state,correlator,noise=1):
    return sym.nsimplify(sym.trace(state*correlator),
                tolerance=1e-10,rational=False).evalf().subs(ALPHA,noise)
    
#%%

#Retrieve the dual values, which will be a dictionary
    
dual_vals = solution.duals

# we are interested in the values corresponding to the probability constraints

inequality_coeff=[]
i=0
for constraint in prob_constraints:
    inequality_coeff.append([dual_vals[constraint][0],prob_coordinates[i]])
    i=i+1

# print the inequality
    
string_ineq=''
suma=0
for ix in inequality_coeff:
    aux_string='p('+str(ix[1][1][0])+str(ix[1][1][1])+str(ix[1][1][2])+'|'+str(ix[1][0][0])+str(ix[1][0][1])+str(ix[1][0][2])+')'
    var=sym.symbols(aux_string)
    # !!
    # now we just round up the numbers and display them nicely, and rescale
    # the inequality
    suma=suma+sym.Rational(sym.Float(ix[0]/inequality_coeff[0][0],2))*var
print("$\sum_i s_{abc}^{xyz}p(abc|xyz)$=")
print(suma)
    






