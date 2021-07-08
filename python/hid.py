from __future__ import division
import picos
import numpy as np
import sympy as sym
import itertools

#from IPython.display import display, Math, Latex # beautiful output

sym.init_printing()

### BIG CONSTANT VALUES
# %%
NR_PARTIES = 2


OUTPUT_VALUES  = [0,1]
NR_OUTPUTS_PER_PARTY = len(OUTPUT_VALUES)

SETTING_VALUES = [0,1]
NR_SETTINGS_PER_PARTY = len(SETTING_VALUES)
###

# %%
'''
In what follows I want to define stuff for defining a quantum state and 
taking measurements on it.
'''

from sympy.physics.quantum import TensorProduct, Dagger

def norm_pure_state(state):
    return sym.sqrt((Dagger(state)*state)[0])

def projector_onto_state(state):
    normm = norm_pure_state(state)
    if normm != 0:    
        return state*Dagger(state)/normm
    else:
        print("ERROR: norm 0")

sig0_norm   = sym.Matrix([[1, 0],[0, 1]])#/sqrt(2)
sig0_norm   = sig0_norm/sym.trace(sig0_norm)


# All these are written in the Z basis
sig1   = sym.Matrix([[0, 1],[1, 0]])#/sqrt(2)
sig2   = sym.Matrix([[0,-sym.I],[sym.I, 0]])#/sqrt(2)
sig3   = sym.Matrix([[1, 0],[0,-1]])#/sqrt(2)
u1base = sym.Matrix(4,1,([sig0_norm,sig1,sig2,sig3]))


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

rho_werner = ((1-ALPHA)*TensorProduct(sig0_norm,sig0_norm) 
            + ALPHA*(max_ent*Dagger(max_ent)))


# Important function, as it defines what measurement we mean by the 
# different values of a, x, b,y etc
def Meas_proj(a,x):
    if x == 0:
        if a == 0:
            return projector_onto_state(sig3_eigenbasis[1])
        if a == 1:
            return projector_onto_state(sig3_eigenbasis[0])
    if x == 1:
        if a == 0:
            return projector_onto_state(sig2_eigenbasis[1])
        if a == 1:
            return projector_onto_state(sig2_eigenbasis[0])
    
# just takes a state and what each party measures and gives the overall
# probability distribution as a function of alpha, the noise
# 2 inputs inputs: rho, [[a,b,c],[x,y,z]]
# returns: p(a,b,c|x,y,z) = \tr (rho M_a|x tensor M_b|y tensor M_c|z) etc.
def probability(rho, output_settings_list):
    outputs = output_settings_list[0]
    settings = output_settings_list[1]
    kronn = Meas_proj(outputs[0],settings[0])
    for i in range(len(outputs)-1):
        kronn = TensorProduct(kronn,Meas_proj(outputs[i+1],settings[i+1]))
    print(kronn)
    return sym.trace(rho*kronn)



# %%
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


#%%        

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

#%%
    
CONST_deterministic_points = (NR_OUTPUTS_PER_PARTY**NR_PARTIES *
                              NR_SETTINGS_PER_PARTY**NR_PARTIES)

q = picos.RealVariable("q", CONST_deterministic_points)
alpha = picos.RealVariable("alpha")


P = picos.Problem()

P.reset()
P.set_objective("max", alpha)  # alpha=0 -> mixed state, alpha=1 -> pure state
                               # I should get alpha=0.34 (p=0.66)


P.add_constraint(alpha >= 0.0)
P.add_constraint(alpha <= 1.0)


#%%


p_vector_as_dict = {}
p_vector_as_dict_string = {}
for setting_det in itertools.product(SETTING_VALUES, repeat=NR_PARTIES):
    setting_det_list = list(setting_det)
    for output_det in itertools.product(OUTPUT_VALUES, repeat=NR_PARTIES):
        output_det_list = list(output_det)
        idx = translate_tuple_to_index([output_det_list,setting_det_list])
        prob = probability(rho_werner,
                        [setting_det_list, output_det_list])
        p_vector_as_dict[idx] = prob
        p_vector_as_dict_string[dic_index_to_probs_string[idx]] = prob


  
p_vector_as_dict_float = {}
p_vector = []
for i in range(CONST_deterministic_points):
    aux = p_vector_as_dict[i]
    p_vector.append(aux)
    p_vector_as_dict_float[i] = aux

# turns out PICOS and Sympy kind of hate each other, I cannot substitute
# the symbol ALPHA in the expressions with a PICOS affive variable, alpha.
# to solve this I will use lambidfy to convert sympy symbolic expressions
# into functions. as a FUNCTION I can input an affine expression and the
# output will be PICOS compatible
p = []
for i in range(CONST_deterministic_points):
    p.append(sym.lambdify(ALPHA,p_vector[i])(alpha))
    
#%%

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
    one a,b which gives 1, and the rest 0. this is what makes it deterministic
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

#%%


for i in range(CONST_deterministic_points):
    # i is component of det behaviours
    # we impose compoennt by component that the convex combination of det
    # behaviorus shoudl equal our fixed probabilitiy distirbution
    suma = 0
    for j in range(CONST_deterministic_points):
        suma = suma + q[j] * list_of_det_behaviours[j][i]
    P.add_constraint(suma == p[i])

positivity_constraint = P.add_constraint(q >= 0.)

normalization_hidden_var = P.add_constraint(q.sum == 1.)

print(P)

#%%

#solution = P.solve(primals=False)
solution = P.solve(solver="mosek")

print("Problem status:",solution.problemStatus)
print("Variable values:")
print("alpha =",alpha.value)
print("q =",([(q[i].value) for i in range(CONST_deterministic_points)]))
