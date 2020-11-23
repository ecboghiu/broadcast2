load('059value.mat');

p_entangled = ProbMultidimArrayInstrumental(NoisyWernerState(0), best_povm, best_channels);
p_uniform = ProbMultidimArrayInstrumental(NoisyWernerState(1), best_povm, best_channels);
[newalpha, bellcoeffs, LPstatus, dual_alpha] = BroadcastInstrumentLP(p_entangled, p_uniform, ins, outs);

bell_in_prob = dispBellCoeffsINSTR(bellcoeffs, ins, outs);

% clean small coefficients and normalize

[C,T]=coeffs(bell_in_prob);
newC = C/max(abs(C));
newC(abs(newC)<1e-6)=0;
newbell=dot(newC,T);

newbell=vpa(newbell,3);
