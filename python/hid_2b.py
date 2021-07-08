# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 14:24:20 2020

@author: embog
"""

import picos
import numpy as np
import itertools


NR_PARTIES = 2
OUTPUT_VALUES  = [[0,1],[0,1,2,3]]
INPUT_VALUES = [[0,1],[0,1,2,3]]

def test_party_nr():
    assert len(OUTPUT_VALUES) == NR_PARTIES

NR_DET_VERTICES = 1
for i in range(NR_PARTIES):
    extra = len(OUTPUT_VALUES[i])**len(INPUT_VALUES[i])
    NR_DET_VERTICES = NR_DET_VERTICES * extra

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

detstrats = detstrats_ab(2,4,2,4)

P = picos.Problem()

q = picos.RealVariable("q",NR_DET_VERTICES)

P.add_constraint(q >= 0)
P.add_constraint(q.sum == 1)


na, nb, nx, ny = 2, 4, 2, 4
prob=[[[[0 for i in range(ny)]for j in range(nx)]for k in range(nb)] for l in range(na)]
for x,y in itertools.product(*INPUT_VALUES):
    for a,b in itertools.product(*OUTPUT_VALUES):
        summ = 0
        for lam in range(NR_DET_VERTICES):
            summ = summ + detstrats[lam][x,y,a,b]*q[lam]
        prob[a][b][x][y] = summ

# marginalized over a
prob_b=[[[0 for i in range(ny)]for j in range(nx)]for k in range(nb)]
for x,y in itertools.product(*INPUT_VALUES):
    for b in OUTPUT_VALUES[1]:
        summ = 0
        for a in OUTPUT_VALUES[0]:
            summ = summ + prob[a][b][x][y]
        prob_b[b][x][y] = summ
        
# We don't need to add non-signaling conditions between A and B because
# this is built into the local model. We do need though to add NS conditions
# between B1 and B2.

#Now we want to translate between (b1,b2) and b:
def to_idx(b1,b2):
    return b2+b1*2

prob_b1=[[[[0 for l in range(2)] for i in range(2)]for j in range(nx)]for k in range(2)]
for x,y1,y2 in itertools.product([0,1],[0,1],[0,1]):
    for b1 in [0,1]:
        summ = 0
        for b2 in [0,1]:
            summ = summ + prob_b[to_idx(b1,b2)][x][to_idx(y1,y2)]
        prob_b1[b1][x][y1][y2] = summ
        
prob_b2=[[[[0 for l in range(2)] for i in range(2)]for j in range(nx)]for k in range(2)]
for x,y1,y2 in itertools.product([0,1],[0,1],[0,1]):
    for b2 in [0,1]:
        summ = 0
        for b1 in [0,1]:
            summ = summ + prob_b[to_idx(b1,b2)][x][to_idx(y1,y2)]
        prob_b2[b2][x][y1][y2] = summ

# NS constraints when marginalizing over Bob1. We impose that the remaining
# distribution of Bob2 should not depent on the input y1 of bob1. We do not
# care about Alice's input x because by using a local distribution we 
# already know that it is non-signalling.
for b2,y2 in itertools.product([0,1],[0,1]):
    idx = 0
    for y1_1, y1_2 in itertools.combinations([0,1],2):
        #combinations('ABCD', 2) --> AB AC AD BC BD CD
        #we use combinations becausee for no signalling we want for every pair 
        #of different y1's to show that they are different, but y1_1 = y1_2
        #is the same as y1_2 = y1_1, so we use combinations instead of iter.product
        #to avoid adding the same equality twice, or redundant equalities
        #such as y1_1 = y1_1.
        P.add_constraint(prob_b2[b2][x][y1_1][y2] ==
                         prob_b2[b2][x][y1_2][y2])
        
#Same for Bob1 
for b1,y1 in itertools.product([0,1],[0,1]):
    idx = 0
    for y2_1, y2_2 in itertools.combinations([0,1],2):
        P.add_constraint(prob_b1[b1][x][y1][y2_1] ==
                         prob_b1[b1][x][y1][y2_2])

P.set_objective(None)


solution = P.solve(solver="mosek")
print("Problem status:",solution.problemStatus)
print()

