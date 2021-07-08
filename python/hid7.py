# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:04:28 2020

@author: embog
"""

# TODO CHANGE FROM np.matrix to np.array

import picos
import numpy as np
import itertools
import gc
from tqdm import tqdm
import copy

import sympy as sym

GLOBAL_PARTIES = [0,1,2] # 0--> A, 1-> B, 2-->C
OUTPUT_VALUES  = [[0,1],[0,1],[0,1]]
INPUT_VALUES = [[0,1,2],[0,1],[0,1]]
NR_PARTIES = len(OUTPUT_VALUES)

def test_party_nr():
	assert len(INPUT_VALUES) == NR_PARTIES
	
#%%
	
'''
Quantum distribution
'''


def projector_onto_state(state):
	normm2=np.linalg.norm(state)
	if normm2 > 1e-10:    
		return np.dot(state,state.conj().T)/normm2 # |psi><psi|/<psi|psi>
	else:
		print("ERROR: norm 0")
		return 0
	
sig0   = np.array([[1, 0],[0, 1]])
# All these are written in the Z basis, |0>, |1>
sig1   = np.array([[0, 1],[1, 0]]) 
sig2   = np.array([[0,-1j],[1j, 0]])
sig3   = np.array([[1, 0],[0,-1]])

e1 = np.array([[1],[0]])
e2 = np.array([[0],[1]])

e1e1 = np.kron(e1,e1)
e1e2 = np.kron(e1,e2)
e2e1 = np.kron(e2,e1)
e2e2 = np.kron(e2,e2)

Phi_minus = 1/np.sqrt(2) * (e1e1 - e2e2)
Psi_plus  = 1/np.sqrt(2) * (e1e2 + e2e1)
Phi_plus  = 1/np.sqrt(2) * (e1e1 + e2e2)


def ret_proj_eigs(obs):
	obs_val, obs_vec = np.linalg.eigh(obs)
	
	ret = []
	for i in range(np.shape(obs_vec)[0]):
		# the following so we have a column vector
		aux_vec = obs_vec[:,i].reshape((len(obs_vec),1)) 
		ret.append(projector_onto_state(aux_vec))
	return ret

def give_party_projectors(operators):
	party_proj = []
	for x in range(len(operators)):
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
		out_pproj.append(give_party_projectors(operator_choice[party]))
	return out_pproj

def initialize_P_proj():
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
	return P_proj

def ini_state_noise(visibility):
	Id4    = np.kron(sig0,sig0)
	ini_state = np.dot(Phi_plus, Phi_plus.conj().T)
	return visibility * ini_state + (1-visibility) * (1.0/4*Id4)

def state_noise(alpha):
	aux_state = ini_state_noise(alpha)
	alp  = np.pi/8
	psi0 =  np.sin(alp) * Phi_minus + np.cos(alp) * Psi_plus
	psi1 = -np.cos(alp) * Phi_minus + np.sin(alp) * Psi_plus
	
	U      = np.dot(psi0 , e1.conj().T) + np.dot(psi1 , e2.conj().T)
	Big_U  = np.kron(sig0,U)
	
	return np.dot(np.dot(Big_U,aux_state),Big_U.conj().T)

def state_noise_PICOS(alpha):
	aux_state = ini_state_noise(alpha)
	alp  = np.pi/8
	psi0 =  np.sin(alp) * Phi_minus + np.cos(alp) * Psi_plus
	psi1 = -np.cos(alp) * Phi_minus + np.sin(alp) * Psi_plus
	
	U      = np.dot(psi0 , e1.conj().T) + np.dot(psi1 , e2.conj().T)
	Big_U  = np.kron(sig0,U)
	Big_U_T = Big_U.conj().T
	
	
	Big_U  = picos.Constant(Big_U)
	Big_U_T = picos.Constant(Big_U_T)
	
	return Big_U * aux_state * Big_U_T


def probability_vis(visibility,inputs, outputs,P_proj):
	kronn = 1
	for i in range(0,NR_PARTIES):
		kronn = np.kron(kronn,P_proj[i][inputs[i]][outputs[i]])
	#print(kronn)
	return np.trace(np.dot(state_noise(visibility),kronn))
	#return sym.nsimplify(ret,tolerance=1e-10,rational=False).evalf()
	
def probability_vis_PICOS(visibility,inputs,outputs,P_proj):
	kronn = P_proj[0][inputs[0]][outputs[0]]
	for i in range(1,NR_PARTIES):
		kronn = picos.kron(kronn,P_proj[i][inputs[i]][outputs[i]])
	return picos.trace(state_noise_PICOS(visibility)*kronn)

def probability_symbolic(ALPHA,inputs,outputs):
	prob = sym.Matrix(probability_vis(ALPHA,inputs,outputs))
	return sym.nsimplify(prob,tolerance=1e-10,rational=False).evalf()

#%%
	
def prepare_for_picos(picos_variable_alpha,operator_choice):
	#ALPHA=sym.symbols('alpha')
	
	temp=[[[[[[0 for c in range(2)] for b in range(2)] for a in range(2)]
				 for z in range(2)] for y in range(2)] for x in range(3)]
	#print("prepare_for_picos: Be patient, this step takes a while.")
	for x_list in tqdm(itertools.product(*INPUT_VALUES)):
		for a_list in itertools.product(*OUTPUT_VALUES):
			x,y,z = x_list[0], x_list[1], x_list[2]
			a,b,c = a_list[0], a_list[1], a_list[2]
			# in the followint, .tolist() because we have a 
			# ndarray(..,dtype=object)
			
			#prob=probability_vis(ALPHA,[x,y,z],[a,b,c],operator_choice)
			prob=probability_vis_PICOS(picos_variable_alpha,[x,y,z],[a,b,c],operator_choice)
			
			# the following to remove small coeffs like 1.04e-17 and so on
			
			
			#prob=sym.simplify(prob)
			#prob=sym.nsimplify(prob,tolerance=1e-14,rational=False).evalf()
			#prob=sym.lambdify(ALPHA,prob)(picos_variable_alpha)
			
			
			temp[x][y][z][a][b][c]=prob
			
	return temp
		
#%%

def get_critical_visibility(initial_state,operator_choice,
							give_bell_coefficients=False):
	# define the deterministic strategies over Alice
	nr_det_points = 8
	det = [[[0 for a in range(len(OUTPUT_VALUES[0]))]
				 for x in range(len(INPUT_VALUES[0]))]
					for i in range(nr_det_points)]
	counter = 0
	for a_tuple in itertools.product(OUTPUT_VALUES[0],
									 repeat=len(INPUT_VALUES[0])):
		for x in INPUT_VALUES[0]:
			det[counter][x][a_tuple[x]] = 1
		counter = counter + 1


	P = picos.Problem()
	alpha = picos.RealVariable("alpha")
	q = picos.RealVariable("q",nr_det_points*2*2*2*2)
	# TODO change this to an index function i---> lam y z b c, easier
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

	# CONSTRAINTS
	positivity_constraints = []
	positivity_constraints = [*positivity_constraints, alpha >= 0,
													   alpha <= 1]
	positivity_constraints = [*positivity_constraints, q>=0]
	P.add_list_of_constraints(positivity_constraints)

	prob_for_picos = prepare_for_picos(alpha,P_proj)
	
	#print("\n prob_for_picos:\n",prob_for_picos)

			   
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
	#print("\n non_signalling_constraints:\n",non_signalling_constraints)
	P.add_list_of_constraints(non_signalling_constraints)

	prob_constraints, prob_coordinates = [], []
	for a,b,c in itertools.product(*OUTPUT_VALUES):
		for x,y,z in itertools.product(*INPUT_VALUES):
			suma = 0
			for i in range(nr_det_points):
				suma = suma + det[i][x][a] * q_vars[i][y][z][b][c]
			
			#prob1 = probability(final_state,[x,y,z],[a,b,c])
			prob1=prob_for_picos[x][y][z][a][b][c]
			
			prob_coordinates.append([[x,y,z],[a,b,c]])
					
			# change symbolic ALPHA to PICOS variable alpha
			#prob2 = (sym.lambdify(ALPHA,prob1))(alpha)
	
			prob_constraints.append(suma == prob1)
	
	#print("\n prob_constraints:\n",prob_constraints)
	P.add_list_of_constraints(prob_constraints)

	# Now solve the problem.
	P.set_objective("max",alpha)
	#P.set_objective(None)
	
	#print("\nPICOS problem:",P)
	solution = P.solve(solver="cvxopt")
	
	#print("\nProblem status:", solution.problemStatus)
	#print("\nNoise tolerance: alpha=", alpha.value)
	
	return_list=[alpha.value]
	
	if give_bell_coefficients==True:
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
			
		#string_ineq=''
		suma=0
		suma_numerica=0
		#aux_coeff = [inequality_coeff[i][0] for i in range(len(inequality_coeff))]
		normalizing_nr = 1.15470053837925#1#min(aux_coeff)
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
			suma=suma+sym.Rational(sym.Float(np.real(ix[0]/normalizing_nr),2))*var
			#suma=suma+ix[0]/normalizing_nr*var
			suma_numerica=suma_numerica+ix[0]/normalizing_nr*probability_vis(alpha.value,[x,y,z],[a,b,c],operator_choice)
			#suma_numerica=suma_numerica+(ix[0]/normalizing_nr*
			#    probability(final_state,[x,y,z],[a,b,c],opertor_choise).subs(ALPHA,alpha_value))
			s_coeff[x][y][z][a][b][c] = ix[0]
			
		suma=sym.nsimplify(suma,tolerance=1e-10,rational=False).evalf()
		#print(r"$\sum s_{abc}^{xyz}p(abc|xyz)$=", suma,"=",suma_numerica)
		return_list.append(np.array(s_coeff))
		return_list.append([suma,suma_numerica])
	return return_list

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
	given_party_outputs = OUTPUT_VALUES[given_party]
	given_party_inputs  = INPUT_VALUES[given_party]
	
	#print("given_party_outputs=",given_party_outputs)
	#print("given_party_inputs=",given_party_inputs)
	
	variable_list = [[0 for i in range(len(given_party_outputs))]
						for j in range(len(given_party_inputs))]
	for a in given_party_outputs:
		for x in given_party_inputs:
			str_name = "(Pi_"+str(given_party)+")^"+str(x)+"_"+str(a)
			#print(str_name)
			variable_list[x][a] = picos.HermitianVariable(str_name,
						 (len(given_party_outputs),len(given_party_outputs)))
			#problem.add_variable(variable_list[x][a])
			#print(variable_list[x][a])

	#%%

	# normalization constriants, sum over outputs given x should be identity
	normalization_constraints = []
	for x in given_party_inputs:
		summ = 0
		for a in given_party_outputs:
			summ = summ + variable_list[x][a]
		aux_dimension = len(given_party_outputs)
		normalization_constraints.append(summ == np.eye(aux_dimension))
	problem.add_list_of_constraints(normalization_constraints)
	
	#print("normalization_constraints:\n", normalization_constraints)

	#%%

	#every variable is a projector, so I'm 
	sdp_constraints = []
	for x in given_party_inputs:
		for a in given_party_outputs:
			sdp_constraints.append(variable_list[x][a] >> 0)
	problem.add_list_of_constraints(sdp_constraints)
	#print("sdp_constraints:\n", sdp_constraints)

	#%%
	# now we construct the objective function
	# here there's the question of whether I should first sum over everything
	# and then take the trace, or sum the traces of each element
	# I think that it should be simpler to create sum over everything first

	state=state_noise(alpha_value)
	state=np.matrix(state).astype(np.complex128)
	state=picos.Constant(state)
	P_proj=initialize_P_proj()
	P_proj_for_PICOS=copy.deepcopy(P_proj) # we ccopy the list
 
	#TODO:PROBABLY I CAN REMOVE THIS STEP
	# Preprocess P_proj
	for party in GLOBAL_PARTIES:
		for x in INPUT_VALUES[party]:
			for a in OUTPUT_VALUES[party]:
				aux_el = P_proj[party][x][a]
				P_proj_for_PICOS[party][x][a] = picos.Constant(np.matrix(
											aux_el
											).astype(np.complex128))

	summ = 0
	for out_list in tqdm(itertools.product(*OUTPUT_VALUES)):
		for in_list in itertools.product(*INPUT_VALUES):
			x, y, z = in_list[0], in_list[1], in_list[2]
			a, b, c = out_list[0], out_list[1], out_list[2]
			coeffbell=s_bell[x][y][z][a][b][c]
			if abs(coeffbell) > 1e-10:
				#print("in_list=",in_list)
				#special boundary case:
				if given_party == 0:
					tensor_product = variable_list[x][a]
				else:
					#sanitized=np.matrix(P_proj_list[0][in_list[0]][out_list[0]])
					#sanitized=sanitized.astype(np.complex128)
					#tensor_product = sanitized
					tensor_product=P_proj_list[given_party][x][a]
					
				for p_idx in range(1,len(GLOBAL_PARTIES)):
					x_idx, a_idx = in_list[p_idx], out_list[p_idx]
					if p_idx != given_party:
						#sanitized=np.matrix(P_proj_list[p_idx][x_idx][a_idx]).astype(np.complex128)
						tensor_product = picos.kron(tensor_product,P_proj_list[p_idx][x_idx][a_idx])
					else:
						tensor_product = picos.kron(tensor_product,variable_list[x_idx][a_idx])
				summ = summ + coeffbell*tensor_product
				
	obj = state*summ 

	obj = (picos.trace(obj)).real      

	problem.set_objective("max", obj)

	sol=problem.solve(solver="mosek", verbosity=2)
	#print(sol.problemStatus)
	objective_value=sol.value
	
	'''
	solution = [[0 for i in range(len(given_party_outputs))]
				   for j in range(len(given_party_inputs))]
	for a in given_party_outputs:
		for x in given_party_inputs:
			solution[x][a] = np.array(variable_list[x][a].value)
	'''

	new_P_proj = copy.deepcopy(P_proj_list)
	for x in given_party_inputs:
		for a in given_party_outputs:
			new_P_proj[given_party][x][a] = np.array(variable_list[x][a].value)
	

	return objective_value, new_P_proj 

def see_saw_over_measurements(input_state, s_bell_coeffs, ini_P_proj):
    party_nr = len(GLOBAL_PARTIES)
    input_P_proj = copy.deepcopy(ini_P_proj)
    print("\n")
    for i in range(10):
        party = i % party_nr
        print(party)
        obj_val, new_P_proj = optimize_over_measurements(input_state, 
                                                         party,
                                                         s_bell_coeffs,
                                                         ini_P_proj)
        input_P_proj = new_P_proj
        print("obj_val:\n",obj_val)
    

if __name__ == "__main__":
	P_proj=initialize_P_proj()
	
	out = get_critical_visibility(None, P_proj, give_bell_coefficients=True)
	alpha_value = out[0]
	s_coeff     = out[1]
	s_symbollic = out[2][0]

	print("\n alpha_value: \n",alpha_value)
	#print("\n s_coeff: \n",s_coeff)
	#print("\n s_symbollic: \n",s_symbollic)
    
	obj_value, new_P_proj = optimize_over_measurements(state_noise_PICOS(alpha_value),
											given_party=0,
											s_bell=np.array(s_coeff),
											P_proj_list=P_proj)
	
	print("old\n",P_proj)
	print("new\n",new_P_proj)
	#see_saw_over_measurements(state_noise_PICOS(alpha_value),s_coeff,P_proj)

