# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:04:28 2020

@author: embog
"""


import picos
import numpy as np
import itertools

import sympy as sym

#from sympy.physics.quantum import Dagger, TensorProduct

GLOBAL_PARTIES = [0,1,2] # 0--> A, 1-> B, 2-->C
GLOBAL_DIC_NUM_TO_PARTY = {0: 'A', 1: 'B', 2: 'C'}
GLOBAL_LIST_NUM_TO_PARTY = ['A','B','C'] # maybe this is faster than the dic??
GLOBAL_DIC_PARTY_TO_NUM = {'A': 0, 'B': 1, 'C':2}


OUTPUT_VALUES  = [[0,1],[0,1],[0,1]]
INPUT_VALUES = [[0,1,2],[0,1],[0,1]]
NR_PARTIES = len(OUTPUT_VALUES)

def test_party_nr():
    assert len(INPUT_VALUES) == NR_PARTIES

    
#%%
    
'''
Quantum distribution
'''

def norm_pure_state(state):
    # [0] in the output because it outputs a 1x1 matrix annoyingly
    return np.sqrt(state.H*state)

def projector_onto_state(state):
    normm2=state.H*state
    if normm2 != 0:    
        return state*state.H/normm2 # |psi><psi|/<psi|psi>
    else:
        print("ERROR: norm 0")
        return 0
    
sig0   = np.matrix([[1, 0],[0, 1]])
# All these are written in the Z basis, |0>, |1>
sig1   = np.matrix([[0, 1],[1, 0]]) 
sig2   = np.matrix([[0,-1j],[1j, 0]])
sig3   = np.matrix([[1, 0],[0,-1]])

e1 = np.matrix([[1],[0]])
e2 = np.matrix([[0],[1]])

e1e1 = np.kron(e1,e1)
e1e2 = np.kron(e1,e2)
e2e1 = np.kron(e2,e1)
e2e2 = np.kron(e2,e2)

Phi_minus = 1/np.sqrt(2) * (e1e1 - e2e2)
Psi_plus  = 1/np.sqrt(2) * (e1e2 + e2e1)
Phi_plus  = 1/np.sqrt(2) * (e1e1 + e2e2)

ini_state = Phi_plus * Phi_plus.H

alp = np.pi/8
psi0 =  np.sin(alp) * Phi_minus + np.cos(alp) * Psi_plus
psi1 = -np.cos(alp) * Phi_minus + np.sin(alp) * Psi_plus

GLOB_U = psi0 * e1.H+ psi1 * e2.H
GLOB_Big_U = np.kron(sig0,GLOB_U)
GLOB_Id4 = np.kron(sig0,sig0)

def ret_proj_eigs(obs):
    obs_val, obs_vec = np.linalg.eigh(obs)
    ret = []
    for i in range(np.shape(obs_vec)[0]):
        ret.append(projector_onto_state(obs_vec[:,i]))
    return ret

def give_party_projectors(party,operators):
    nr_inputs = len(INPUT_VALUES[party])
    if len(operators) != nr_inputs:
        exit('Not enough operators')
    party_proj = []
    for x in range(nr_inputs):
        party_proj.append(ret_proj_eigs(operators[x]))
    return party_proj

def generate_P_proj(operator_choice):
    '''
    Auxiliary function. It takes as input "operator_choice" which is 
    something of the form [[A0,A1,A2],[B0,B1],[C0,C1]] ie, operators for 
    measurements for each party and it outputs the same list but furter
    decomposing the operators into projectors onto its eigenspaces, like
    A0 --> [Proj_eig1_A0, proj_eig2_A0] etc. 
    '''
    out_pproj=[]
    for party in GLOBAL_PARTIES:
        out_pproj.append(give_party_projectors(party,operator_choice[party]))
    return out_pproj

#Alice
A0 = sig3 #sig_z
A1 = sig1 #sig_x
A2 = sig2 #sig_y

#Bob
phi = np.arctan(1/np.sqrt(2))
B0 = np.cos(phi) * sig1 + np.sin(phi) * sig2
B1 = np.cos(phi) * sig1 - np.sin(phi) * sig2

#Charie
C0 = sig3
C1 = sig1

P_proj = generate_P_proj([[A0,A1,A2],[B0,B1],[C0,C1]])


#%%

#U

def ini_state_noise(visibility):
    return visibility * ini_state + (1-visibility) * (1.0/4*GLOB_Id4)

def state_noise(alpha):
    aux_state = ini_state_noise(alpha)
    return GLOB_Big_U * aux_state * GLOB_Big_U.H

def meas_proj(party,i,o):
    return P_proj[party][i][o]

def probability_vis(visibility,inputs, outputs):
    kronn = 1
    for i in range(0,NR_PARTIES):
        kronn = np.kron(kronn,meas_proj(i,inputs[i],outputs[i]))
    #print(kronn)
    ret = np.real(np.trace(state_noise(visibility)*kronn))
    return ret
    #return sym.nsimplify(ret,tolerance=1e-10,rational=False).evalf()

def probability_symbolic(ALPHA,inputs,outputs):
    prob = sym.Matrix(probability_vis(ALPHA,inputs,outputs))
    return sym.nsimplify(prob,tolerance=1e-10,rational=False).evalf()
        
#%%

nr_det_points = 8
det = [[[0 for a in range(len(OUTPUT_VALUES[0]))]
             for x in range(len(INPUT_VALUES[0]))]
                for i in range(nr_det_points)]

counter = 0
for a_tuple in itertools.product(OUTPUT_VALUES[0],
                                 repeat=len(INPUT_VALUES[0])):
    a_list = list(a_tuple)
    for x in INPUT_VALUES[0]:
        det[counter][x][a_tuple[x]] = 1
    counter = counter + 1


P = picos.Problem()


#%%
# CONSTRAINTS
positivity_constraints = []


alpha = picos.RealVariable("alpha")
positivity_constraints = [*positivity_constraints, alpha >= 0,
                                                   alpha <= 1]


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

final_state = state_noise(alpha)

for a,b,c in itertools.product(*OUTPUT_VALUES):
    for x,y,z in itertools.product(*INPUT_VALUES):
        
        suma = 0
        for i in range(nr_det_points):
            suma = suma + det[i][x][a] * q_vars[i][y][z][b][c]
        
        prob1 = probability(final_state,[x,y,z],[a,b,c])
        prob_coordinates.append([[x,y,z],[a,b,c]])
                
        # change symbolic ALPHA to PICOS variable alpha
        #prob2 = (sym.lambdify(ALPHA,prob1))(alpha)

        prob_constraints.append(suma == prob1)

P.add_list_of_constraints(prob_constraints)


#%%
# Now solve the problem.

P.set_objective("max",alpha)
#P.set_objective(None)

solution = P.solve(solver="mosek")
print("Problem status:", solution.problemStatus)

print("Noise tolerance: alpha=", alpha.value)

solution
#%%

'''
This is to compare with the things from Joe's notes.
'''

'''
correlator = ( -2*TensorProduct(A2,B0,sig0)+2*TensorProduct(A2,B1,sig0)
              + TensorProduct(A0,B0,C0) + TensorProduct(A0,B1,C1)
              + TensorProduct(A1,B1,C1) - TensorProduct(A1,B0,C0)
              + TensorProduct(A0,B0,C1) + TensorProduct(A0,B1,C0)
              + TensorProduct(A1,B0,C1) - TensorProduct(A1,B1,C0) )

def corr_noise_resistance(state,correlator,noise=1):
    return np.trace(state,correlator)
    #return sym.nsimplify(sym.trace(state*correlator),
     #           tolerance=1e-10,rational=False).evalf().subs(ALPHA,noise)
  '''  
#%%

#Retrieve the dual values, which will be a dictionary
    
dual_vals = solution.duals

# we are interested in the values corresponding to the prob constraints


s_coeff = [[[[[[0 for c in range(len(OUTPUT_VALUES[2]))]
                  for b in range(len(OUTPUT_VALUES[1]))]
                  for a in range(len(OUTPUT_VALUES[0]))]
                  for z in range(len(INPUT_VALUES[2]))]
                  for y in range(len(INPUT_VALUES[1]))]
                  for x in range(len(INPUT_VALUES[0]))]

inequality_coeff=[]
i=0
for constraint in prob_constraints:
    inequality_coeff.append([dual_vals[constraint][0],prob_coordinates[i]])
    i=i+1

# print the inequality
    
string_ineq=''
suma=0
for ix in inequality_coeff:
    x,y,z = ix[1][0][0], ix[1][0][1], ix[1][0][2]
    a,b,c = ix[1][1][0], ix[1][1][1], ix[1][1][2]
    s_coeff[x][y][z][a][b][c] = ix[0]
    aux_string=('p('+str(ix[1][1][0])+str(ix[1][1][1])+str(ix[1][1][2])+'|'+
                     str(ix[1][0][0])+str(ix[1][0][1])+str(ix[1][0][2])+')')
    var=sym.symbols(aux_string)
    # !!
    # now we just round up the numbers and display them nicely, and rescale
    # the inequality
    suma=suma+sym.Rational(sym.Float(np.real(ix[0]/
                                             inequality_coeff[0][0]),2))*var
    s_coeff[x][y][z][a][b][c] = sym.Rational(sym.Float(np.real(ix[0]/
                                             inequality_coeff[0][0]),2))
    
print("$\sum s_{abc}^{xyz}p(abc|xyz)$=", suma)

#%%
'''
# the complete expression involves:
# the quantum state with a fixed noise + 
# choices of measurement operators for Alice, Bob, Charlie + 
# a specific quantum channel
# we need to optimize over one while keeping the other ones fixed

#def max_sdp(input_H,dims):
input_H=np.kron(np.kron(sig2,sig1),sig3)
problem = picos.Problem()
H = picos.Constant(input_H)
P = picos.HermitianVariable("P", H.shape)   # if P is PSD then it must 
                                            # be also Hermitian but not the
                                            # other way around
problem.set_objective("max", picos.trace(H*P))
problem.add_constraint(P >> 0) 
problem.add_constraint(P.partial_trace(0,[4,2]) == np.eye(2))
problem.solve(solver="mosek")

print('made it here')
'''
#%%

#Optimize over measurements.

# assume a fixed value of alpha

def optimize_over_measurements(state,given_party,s_bell,P_proj_list):
    r'''
    This runs an optimization by fixing all measurements of the Bell
    expression given by s_bell for all parties except given_party. This gives
    something of the form
    
    .. math::
        \sum_{xyzabc}\operatorname{tr}\rho s^{xyz}_{abc}\Pi^x_a\Pi^y_b\Pi^z_c
        
    where if say given_party==1, then since 1 is 'B', we are taking as 
    SDP variables $\Pi^y_b$ and taking \Pi^x_a,\Pi^z_c as fixed. Then if take
    as an SDP objective the above expression, it is linear and it can be 
    solved. We add as constraints that \sum_b \Pi^y_b = Id \forall y. The
    values for \Pi^x_a,\Pi^z_c in this example are given in P_proj_list. The
    coefficients s^{xyz}_{abc} are given by >>>s_bell[x][y][z][a][b][c]. The 
    quantum state \rho is given by `state`.
    
    state: numpy positive semi definite matrix
        The quantum state, which is simply a Postive-Semidefinite operatos
        with trace equal to one.
    
    given_party: int 
        Tells us the measurement projects of which party we will take as 
        variable and the rest constant. This means that we will ignore 
        in P_proj_list the operatosr corresponding to "party".
            
    s_bell: 6-levels nested list
        A nested list such that s[x][y][z][a][b][c] gives the expression
        of the coefficient of p(abc|xyz) in the Bell inequality.
            
    P_proj_list: 3 levels nested list
        A nested list such that P_proj_list[party][input][output]
        gives the projector onto the specified output conditioned on the 
        input on the party we are interested in.
    '''
    problem = picos.Problem()
    
    #  First define the variables we will use.
    given_party_outputs=OUTPUT_VALUES[given_party]
    given_party_inputs=INPUT_VALUES[given_party]
    print("given_party_outputs=",given_party_outputs)
    print("given_party_inputs=",given_party_inputs)
    variable_list = [[0 for i in range(len(given_party_outputs))]
                        for j in range(len(given_party_inputs))]
    for a in given_party_outputs:
        for x in given_party_inputs:
            str_name = "(Pi_"+str(given_party)+")^"+str(x)+"_"+str(a)
            #print(str_name)
            variable_list[x][a] = picos.HermitianVariable(str_name,
                         len(given_party_outputs))
            #problem.add_variable(variable_list[x][a])
            #print(variable_list[x][a])
    
    # normalization constriants, sum over outputs given x should be identity
    normalization_constraints = []
    for x in given_party_inputs:
        summ = 0
        for a in given_party_outputs:
            summ = summ + variable_list[x][a]
        aux_dimension = len(given_party_outputs)
        normalization_constraints.append(summ == np.eye(aux_dimension))
    print("normalization_constraints",normalization_constraints)
    
    # now we construct the objective function
    # here there's the question of whether I should first sum over everything
    # and then take the trace, or sum the traces of each element
    # I think that it should be simpler to create sum over everything first
    summ = 0
    for out_list in itertools.product(*OUTPUT_VALUES):
        print("out_list=",out_list)
        for in_list in itertools.product(*INPUT_VALUES):
            print("in_list=",in_list)
            #special boundary case:
            if given_party == 0:
                tensor_product = variable_list[in_list[0]][out_list[0]]
                for p_idx in range(1,len(GLOBAL_PARTIES)):
                    x = in_list[p_idx]
                    a = out_list[p_idx]
                    tensor_product = picos.kron(tensor_product,
                                                P_proj_list[p_idx][x][a])
            else:
                tensor_product = P_proj_list[0][in_list[0]][out_list[0]]
                for p_idx in range(1,len(GLOBAL_PARTIES)):
                    x = in_list[p_idx]
                    a = out_list[p_idx]
                    if p_idx != given_party:
                        tensor_product = picos.kron(tensor_product,
                                                    P_proj_list[p_idx][x][a])
                    else:
                        tensor_product = picos.kron(tensor_product,
                                                    variable_list[x][a])
            # this should give us the correct
            # \Pi^x_a \otimes \Pi_y_b \otimes \Pi_z_c
            summ = summ + float(s_bell[x][y][z][a][b][c]) * tensor_product
    objective = picos.trace(state*summ)
    
    problem.set_objective("max", objective)
    
    problem.solve(solver="mosek")
    
    solution = [[0 for i in range(len(given_party_outputs))]
                   for j in range(len(given_party_inputs))]
    for a in given_party_outputs:
        for x in given_party_inputs:
            solution[x][a] = variable_list[x][a].value
    
    return solution

#solution = optimize_over_measurements(state_noise(alpha.value), 1, s_coeff, P_proj)