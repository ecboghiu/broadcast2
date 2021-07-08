# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 15:46:30 2020

@author: embog
"""

import picos
import numpy as np

import sympy as sym
import itertools

NR_PARTIES = 2

OUTPUT_VALUES  = [0,1]
NR_OUTPUTS_PER_PARTY = len(OUTPUT_VALUES)

SETTING_VALUES = [0,1]
NR_SETTINGS_PER_PARTY = len(SETTING_VALUES)

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

sig1_basis2, sig1_eigenvals2 = sig1.diagonalize()
sig1_eigenbasis = [ sig1_basis2[:,i]/norm_pure_state(sig1_basis2[:,i])
                     for i in range((sig1_basis2).shape[0])]
sig1_eigenvalues = [ sig1_eigenvals2[i,i]
                     for i in range((sig1_eigenvals2).shape[0])]


sig2_basis2, sig2_eigenvals2 = sig2.diagonalize()
sig2_eigenbasis = [ sig2_basis2[:,i]/norm_pure_state(sig2_basis2[:,i])
                     for i in range((sig2_basis2).shape[0])]
sig2_eigenvalues = [ sig2_eigenvals2[i,i]
                     for i in range((sig2_eigenvals2).shape[0])]


sig3_basis2, sig3_eigenvals2 = sig3.diagonalize()
sig3_eigenbasis = [ sig3_basis2[:,i]/norm_pure_state(sig3_basis2[:,i])
                     for i in range((sig3_basis2).shape[0])]
sig3_eigenvalues = [ sig3_eigenvals2[i,i]
                     for i in range((sig3_eigenvals2).shape[0])]

### IMPORTANT CONSTANT
ALPHA = sym.symbols('alpha')

e1 = sym.Matrix([1,0])
e2 = sym.Matrix([0,1])

max_ent = TensorProduct(e1,e1)+TensorProduct(e2,e2)
max_ent = max_ent / norm_pure_state(max_ent) 

rho_werner = ((1-ALPHA)*TensorProduct(sig0/sym.trace(sig0),sig0/sym.trace(sig0)) 
            + ALPHA*(max_ent*Dagger(max_ent)))



# IMPORTANT: THIS IS WHERE I DEFINE THE "DIRECTIONS" ON THE BLOCH SPHERE FOR MY 
# MEASUREMENT OPERATORS. N1, N2 for Alice, M1, M2 for Bob.
# No need to normalize them.

N1 = np.array([1,0,0])
N2 = np.array([0,1,0])

M1 = np.array([1,1,1])
M2 = np.array([-1,1,0])
'''

np.random.seed(0)

N1 =  -1+2*np.random.rand(3)
N1 = N1/np.sqrt(np.dot(N1,N1))

N2 =  -1+2*np.random.rand(3)
N2 = N2/np.sqrt(np.dot(N2,N2))

M1 =  -1+2*np.random.rand(3)
M1 = M1/np.sqrt(np.dot(M1,M1))

M2 =  -1+2*np.random.rand(3)
M2 = M2/np.sqrt(np.dot(M2,M2))
'''
print('N1',N1, '\nN2',N2, '\nM1',M1, '\nM2',M2)

A_obs1 = N1[0] * sig1 + N1[1] * sig2 + N1[2] * sig3
A_obs1_eigbasis, A_obs1_eigvals = A_obs1.diagonalize()
    
A_obs2 = N2[0] * sig1 + N2[1] * sig2 + N2[2] * sig3
A_obs2_eigbasis, A_obs2_eigvals = A_obs2.diagonalize()

B_obs1 = M1[0] * sig1 + M1[1] * sig2 + M1[2] * sig3
B_obs1_eigbasis, B_obs1_eigvals = B_obs1.diagonalize()
        
B_obs2 = M2[0] * sig1 + M2[1] * sig2 + M2[2] * sig3
B_obs2_eigbasis, B_obs2_eigvals = B_obs2.diagonalize()

# Important function, as it defines what measurement we mean by the 
# different values of a, x, b,y etc
def Meas_projA(a,x):
    if x == 0:
        if a == 0:
            #return projector_onto_state(e1 + sym.exp(sym.I*phi)*e2)
            return projector_onto_state(A_obs1_eigbasis[:,0])
        elif a == 1:
            #return projector_onto_state(e1 - sym.exp(sym.I*phi)*e2)
            return projector_onto_state(A_obs1_eigbasis[:,1])
        else:
            print("error not the case we want 1")
            return 0
    elif x == 1:
        if a == 0:
            #return projector_onto_state(e1 + sym.exp(sym.I*phi)*e2)
            return projector_onto_state(A_obs2_eigbasis[:,0])
        elif a == 1:
            #return projector_onto_state(e1 - sym.exp(sym.I*phi)*e2)
            return projector_onto_state(A_obs2_eigbasis[:,1])
        else:
            print("error no the case we want 2")
            return 0
    else:
        print("error not the case we want 3")
        return 0
        
def Meas_projB(b,y):
    if y == 0:
        if b == 0:
            #return projector_onto_state(e1 + sym.exp(sym.I*phi)*e2)
            return sym.simplify(projector_onto_state(B_obs1_eigbasis[:,0]))
        elif b == 1:
            #return projector_onto_state(e1 - sym.exp(sym.I*phi)*e2)
            return sym.simplify(projector_onto_state(B_obs1_eigbasis[:,1]))
        else:
            print("error not the case we want 1")
            return 0
    elif y == 1:
        if b == 0:
            #return projector_onto_state(e1 + sym.exp(sym.I*phi)*e2)
            return sym.simplify(projector_onto_state(B_obs2_eigbasis[:,0]))
        elif b == 1:
            #return projector_onto_state(e1 - sym.exp(sym.I*phi)*e2)
            return sym.simplify(projector_onto_state(B_obs2_eigbasis[:,1]))
        else:
            print("error no the case we want 2")
            return 0
    else:
        print("error not the case we want 3")
        return 0
    
# just takes a state and what each party measures and gives the overall
# probability distribution as a function of alpha, the noise
# 2 inputs inputs: rho, [[a,b,c],[x,y,z]]
# returns: p(a,b,c|x,y,z) = \tr (rho M_a|x tensor M_b|y tensor M_c|z) etc.
'''
def probability(rho, output_settings_list):
    outputs = output_settings_list[0]
    settings = output_settings_list[1]
    kronn = Meas_proj(outputs[0],settings[0])
    for i in range(len(outputs)-1):
        kronn = TensorProduct(kronn,Meas_proj(outputs[i+1],settings[i+1]))
    print(kronn)
    return sym.trace(rho*kronn)
'''

def probability(rho, output_settings_list):
    outputs = output_settings_list[0]
    settings = output_settings_list[1]
    kronn = TensorProduct(Meas_projA(outputs[0],settings[0]),
                          Meas_projB(outputs[1],settings[1]))
    #return sym.nsimplify(sym.simplify(sym.trace(rho*kronn)),tolerance=1e-10,rational=True).evalf()
    return sym.simplify(sym.trace(rho*kronn).evalf())
    
P = picos.Problem()

q = picos.RealVariable("q",16)
alpha = picos.RealVariable("alpha",1)


def binary_to_idx(a,b,x,y):
    return y + x * 2 + b * (2**2)+ a * (2**3)


pr = []
dic = {}
idx = 0
for a,b,x,y in itertools.product([0,1],repeat=4):
    pr.append((sym.lambdify(ALPHA,probability(rho_werner, [[a,b], [x,y]])))(alpha))
    dic[idx]=[a,b,x,y]
    idx = idx + 1


def list_to_string(lista):
    '''the lista input should be in the form [[inputs],[settings]]'''
    aux_string = ('p('+
                       (''.join(str(lista[0][e])+','
                                     for e in range(len(lista[0])-1)))
                       +str(lista[0][-1]) # took last element from 
                                       #the loop because the extra 
                                       # comma gave me lsits like 
                                       # [0,2,3,] which is annoying
                                       
                       +'|'+(''.join(str(lista[1][e])+','
                                     for e in range(len(lista[1])-1)))
                       +str(lista[1][-1]) # took last element from 
                                       #the loop because the extra 
                                       # comma gave me lsits like 
                                       # [0,2,3,] which is annoying
                       +')')
    return aux_string

# WARNING: 
# WARNING:
# does not work with double digit  integers
def string_to_list(string):
    '''string should be in the format 'p(i,j,k|m,n,o)' 
    outputs a list such as [[i,j,k],[m,n,o]]
    DOES NOT WORK WITH DOUBLE DIGIT INTEGERS
    thats quite complicated numerically so I don't forsee
    getting close to such complicated cases'''

    # note we assume first two strings are 'p('
    # now we loop for ijk until we hit '|' and then the others until we hit
    # ')'
    idx     = 0
    aux_str = string[idx]
    while aux_str != '(':
        idx     = idx + 1
        aux_str = string[idx]
    
    aux_str  = string[idx]
    aux_str2 = string[idx]
    while aux_str != '|':
        idx      = idx + 1
        aux_str  = string[idx]
        aux_str2 = aux_str2 + string[idx]

    lista_aux = list(aux_str2)
    lista_aux.pop( 0)
    lista_aux.pop(-1)
    len_aux = lista_aux.count(',')
    for i in range(len_aux):
        lista_aux.remove(',')
        
    aux_str  = string[idx]
    aux_str2 = string[idx]

    while aux_str != ')':
        idx      = idx + 1
        aux_str  = string[idx]
        aux_str2 = aux_str2 + string[idx]
    lista_aux2 = list(aux_str2)
    lista_aux2.pop(0)
    lista_aux2.pop(-1)
    len_aux = lista_aux2.count(',')
    for i in range(len_aux):
        lista_aux2.remove(',')
    
    inputs   = [int(i) for i in lista_aux ]
    settings = [int(i) for i in lista_aux2]
    
    return [inputs,settings]
    
# %%
    
dic_index_to_probs_list = {}
dic_index_to_probs_string = {}
dic_probs_to_index = {}  # this will have lists as 'p(a,b,c|x,y,z)' as a  
                         # string because dics cannot have keys as lists,
                         # but only as strings

aux_list = []
idx = 0
for output_det in itertools.product(OUTPUT_VALUES, repeat=NR_PARTIES):
    for setting_det in itertools.product(SETTING_VALUES, repeat=NR_PARTIES):
        aux_list.append([list(output_det),list(setting_det)])
        dic_index_to_probs_list[idx] = aux_list[-1]
        aux_string = list_to_string(aux_list[-1])
        dic_index_to_probs_string[idx] = aux_string
        dic_probs_to_index[aux_string] = idx
        idx = idx + 1       

# useful for chaning from [party,outcome,setting] to a single number label
# eg: [0,0,0] --> 0, [0,1,0] -->2, [0,2,0]--> 3, [1,0,0] --> 5 etc.
def translate_tuple_to_index(input_list):
    '''
    Input list will be of the form [[1,2,3],[2,3,4]]=[inputs,settings] if
    there are 3 parties
    output is the index it corresponds with
    '''
    lista = input_list[0] + input_list[1]
    lista_long = len(lista)
    suma = 0
    for i in range(lista_long):
        suma = suma + lista[lista_long-1-i] * NR_OUTPUTS_PER_PARTY ** (i)
    return suma

CONST_deterministic_points = (NR_OUTPUTS_PER_PARTY**NR_PARTIES *
                              NR_SETTINGS_PER_PARTY**NR_PARTIES) # == 16


p_vector_as_dict = {}
p_vector_as_dict_string = {}
for setting_det in itertools.product(SETTING_VALUES, repeat=NR_PARTIES):
    setting_det_list = list(setting_det)
    for output_det in itertools.product(OUTPUT_VALUES, repeat=NR_PARTIES):
        output_det_list = list(output_det)
        idx = translate_tuple_to_index([output_det_list,setting_det_list])
        prob = (probability(rho_werner,
                        [setting_det_list, output_det_list]))
        p_vector_as_dict[idx] = prob
        p_vector_as_dict_string[dic_index_to_probs_string[idx]] = prob

  
p_vector_as_dict_float = {}
p_vector = []
for i in range(CONST_deterministic_points):
    aux = sym.simplify(p_vector_as_dict[i])
    p_vector.append((aux))
    p_vector_as_dict_float[i] = aux

p = []
for i in range(CONST_deterministic_points):
    p.append(sym.lambdify(ALPHA,p_vector[i])(alpha))

def give_det_behaviour(lam):
    '''
    the structures should be lam=[lam1,lam2,...,lamNR_PARTIES]
    with lam1 for eg being [[0,a1],[1,a2],...,[m,am]] which is {[x,a_x]} for x
    settings or "inputs" and the a's being outputs. For each party, for every 
    possible setting we establish an output because this is a completely 
    deterministic behaviour.
    the output is a vector q_lam such that for component corresponding
    to index, say, (ab|xy) the component is 0 except for when
    a == a_x and b == b_y, where a_x, b_y are from lam, [x,a_x] and [y,b_y]
    so lam determines what q(ab|xy) are =1, and for each x,y there is only 
    one a,b which gives 1, and the rest 0.
    '''

    result = np.zeros(CONST_deterministic_points)

    for setting_det in itertools.product(SETTING_VALUES, repeat=NR_PARTIES):
        setting_det_list = list(setting_det)
        output_det_list = []
        idx_party = 0
        for m in setting_det_list:
            party_det_list = lam[idx_party]
            for pair in party_det_list:
                if pair[0] == m:
                    output_det_list.append(pair[1])
                    break
            idx_party = idx_party + 1
        idx = translate_tuple_to_index([output_det_list,setting_det_list])
        result[idx] = 1
        
    return result

def print_prob_dis(prob_vector):
    for i in range(len(prob_vector)):
        print(dic_index_to_probs_string[i]+' = '+str(prob_vector[i]))


list_of_lambdas = []
short_list_of_lambdas = []
for outs_per_setting in itertools.product(itertools.product(OUTPUT_VALUES,
                                            repeat=NR_SETTINGS_PER_PARTY),
                                            repeat=NR_PARTIES):
    aux_all_parties_with_settings = []
    for party in range(NR_PARTIES):
        aux_all_settings_for_a_party = []
        for sett in range(NR_SETTINGS_PER_PARTY):
            aux_all_settings_for_a_party.append([sett,
                                        outs_per_setting[party][sett]])
        aux_all_parties_with_settings.append(aux_all_settings_for_a_party)
    short_list_of_lambdas.append([list(outs_per_setting[0]),list(outs_per_setting[1])])
    list_of_lambdas.append(aux_all_parties_with_settings)


list_of_det_behaviours = []
for i in range(CONST_deterministic_points):
    list_of_det_behaviours.append(give_det_behaviour(list_of_lambdas[i]))
    
list_of_dics_of_det_behaviours = []
list_of_dics_of_det_behaviours_string = []
for i in range(CONST_deterministic_points):
    aux_dic = {}
    aux_dic_string = {}
    for j in range(CONST_deterministic_points):
        auxx = list_of_det_behaviours[i][j]
        aux_dic[j] = auxx
        aux_dic_string[dic_index_to_probs_string[j]] = auxx
    list_of_dics_of_det_behaviours.append(aux_dic)
    list_of_dics_of_det_behaviours_string.append(aux_dic_string)

def detstrats_ab(nx,ny,na,nb):
    ## function plagiarized from joe
    #create all deterministic strategies in the standard Bell scenario
    #output is numpy multidim array: e.g. detstrats[i,:,:,:,:] gives an np array for
    #the i'th deterministic strategy in the form p(x,y,a,b) <- NOTE NOT p(a,b,x,y) [I prefer the output like this]

    siz1=[na for i in range(nx)]
    D1=[[0 for i in range(nx)] for j in range(na**nx)]
    D1=np.array(D1)
    #convert index to a deterministic strategy for Alice
    #D1[i][:] gives the value of a Alice outputs for each of her inputs for the strategy i
    for i in range(na**nx):
        D1[i][:]=np.array(np.unravel_index(i,siz1))
    #Same for Bob with D2
    D2=[[0 for i in range(ny)] for j in range(nb**ny)]
    D2=np.array(D2)
    if na==nb and nx==ny:
        D2=D1
    else:
        siz2=[nb for i in range(ny)]
        for i in range(nb**ny):
            D2[i][:]=np.array(np.unravel_index(i,siz2))
    #Create p(a,b,x,y) object to store 
    dstrats=np.array([[[[[0 for i in range(nb)]for j in range(na)]for k in range(ny)] for l in range(nx)] for m in range(na**nx*nb**ny)])
    #fill probability table
    c=0
    for cA in range(na**nx):
        for cB in range(nb**ny):
            for x in range(nx):
                for y in range(ny):
                    dstrats[c,x,y,D1[cA][x],D2[cB][y]]=1
            c=c+1
    return dstrats

detstrats = detstrats_ab(2,2,2,2)

for i in range(CONST_deterministic_points):
    # i is component of det behaviours
    # we impose compoennt by component that the convex combination of det
    # behaviorus shoudl equal our fixed probabilitiy distirbution
    idx_list = dic_index_to_probs_list[i]
    a = idx_list[0][0]
    b = idx_list[0][1]
    x = idx_list[1][0]
    y = idx_list[1][1]
    suma = 0
    for j in range(CONST_deterministic_points):
        suma = suma + q[j] * list_of_det_behaviours[j][i]
        #suma = suma + q[j] * detstrats[j,x,y,a,b]
    prob = probability(rho_werner, [[a,b], [x,y]])
    prob = (sym.lambdify(ALPHA,prob))(alpha)
    P.add_constraint(suma == prob)

'''
CONST_deterministic_points=16
for setting_det in itertools.product(SETTING_VALUES, repeat=NR_PARTIES):
    setting_det_list = list(setting_det)
    for output_det in itertools.product(OUTPUT_VALUES, repeat=NR_PARTIES):
        output_det_list = list(output_det)
        x = setting_det_list[0]
        y = setting_det_list[1]
        a = output_det_list[0]
        b = output_det_list[1]
        prob = probability(rho_werner, [[a,b], [x,y]])
        suma = 0
        for j in range(CONST_deterministic_points):
            suma = suma + q[j] * detstrats[j,x,y,a,b]
        proba = (sym.lambdify(ALPHA,prob))(alpha)
        P.add_constraint(suma == proba)
'''
        
'''

for a,b in itertools.product([0,1],repeat=2):
    summ = 0
    for x,y in itertools.product([0,1],repeat=2):
        summ = summ + q[binary_to_idx(a,x,b,y)]
    P.add_constraint(summ == pr[binary_to_idx(a,b,0,0)])        

for a,b in itertools.product([0,1],repeat=2):
    summ = 0
    for x,y in itertools.product([0,1],repeat=2):
        summ = summ + q[binary_to_idx(a,x,y,b)]
    P.add_constraint(summ == pr[binary_to_idx(a,b,0,1)])   
    
for a,b in itertools.product([0,1],repeat=2):
    summ = 0
    for x,y in itertools.product([0,1],repeat=2):
        summ = summ + q[binary_to_idx(x,a,b,y)]
    P.add_constraint(summ == pr[binary_to_idx(a,b,1,0)])   
    
for a,b in itertools.product([0,1],repeat=2):
    summ = 0
    for x,y in itertools.product([0,1],repeat=2):
        summ = summ + q[binary_to_idx(x,a,y,b)]
    P.add_constraint(summ == pr[binary_to_idx(a,b,1,1)])   
'''

P.set_objective("max", alpha)

q_mas   = (1+alpha/np.sqrt(2))/4
q_menos = (1-alpha/np.sqrt(2))/4

'''
P.add_constraint( q[0]  + q[1]  + q[5]  + q[4]   == q_menos) #
P.add_constraint( q[4]  + q[5]  + q[13] + q[12]  == q_menos) #
P.add_constraint( q[0]  + q[2]  + q[8]  + q[10]  == q_mas) #
P.add_constraint( q[4]  + q[6]  + q[12] + q[14]  == q_menos) #
P.add_constraint( q[8]  + q[9] + q[13] + q[12]  == q_mas) #
P.add_constraint( q[0]  + q[2]  + q[6]  + q[4]   == q_mas)  #
P.add_constraint( q[8]  + q[10] + q[12] + q[14]  == q_menos) #
P.add_constraint( q[0]  + q[1]  + q[9] + q[8]   == q_mas) #
P.add_constraint( q[2]  + q[3]  + q[7]  + q[6]   == q_mas)  #
P.add_constraint( q[6]  + q[7]  + q[14] + q[15]  == q_mas) #
P.add_constraint( q[1]  + q[3]  + q[9] + q[11]  == q_menos) #
P.add_constraint( q[5]  + q[7]  + q[13] + q[15]  == q_mas) #
P.add_constraint( q[10] + q[11] + q[14] + q[15]  == q_menos)  #
P.add_constraint( q[1]  + q[3]  + q[7]  + q[5]   == q_menos) #
P.add_constraint( q[9] + q[11] + q[13] + q[15]  == q_mas) #
P.add_constraint( q[2]  + q[3]  + q[11] + q[10]  == q_menos) #
'''

P.add_constraint(q >= 0)
P.add_constraint(q <= 1)

P.add_constraint(alpha >= 0)
P.add_constraint(alpha <= 1)


P.add_constraint(q.sum == 1)

print(P)

solution = P.solve(solver="mosek")

print("Problem status:",solution.problemStatus)
print("Variable values:")

print("alpha  =", alpha.value)
print("q =", q.value)