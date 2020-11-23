mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

B1 = sym('B1');
B2 = sym('B2');
A1 = sym('A1');
A2 = sym('A2');

p_11_11 = sym('p_11_11');
p_11_12 = sym('p_11_12');
p_11_21 = sym('p_11_21');
p_11_22 = sym('p_11_22');

p_12_11 = sym('p_12_11');
p_12_12 = sym('p_12_12');
p_12_21 = sym('p_12_21');
p_12_22 = sym('p_12_22');

p_21_11 = sym('p_21_11');
p_21_12 = sym('p_21_12');
p_21_21 = sym('p_21_21');
p_21_22 = sym('p_21_22');

p_22_11 = sym('p_22_11');
p_22_12 = sym('p_22_12');
p_22_21 = sym('p_22_21');
p_22_22 = sym('p_22_22');

FLAG_Use01obsInsteadOfCorrelator = true;

ins = [2,2];
outs = [2,2];
CHSH = A1*B1+A1*B2+A2*B1-A2*B2-A1-B1;
PAPER = CHSH;
PAPERinprob = ToProbabilityNotationIneqSym(PAPER,ins,outs,FLAG_Use01obsInsteadOfCorrelator);
%PAPERinprob = p_11_11 + p_11_12 + p_11_21 - p_11_22 - ( p_11_11+p_12_11 + p_11_11+p_21_11 );
PAPERbellcoeffs = GetBellCoeffsFromProbSymIneq(PAPERinprob,ins,outs);
LBoundPAPER_NUM = 0;
localbound = ClassicalOptInequalityNoBroadcast(PAPERbellcoeffs)

Pproj = givePprojRANDmaxCHSH();

e1 = [1;0];
e2 = [0;1];
e1e2 = Tensor(e1,e2);
e2e1 = Tensor(e2,e1);
Phi_minus = 1/2^0.5 * (e1e2 - e2e1);
CHSHstate = Phi_minus * Phi_minus';

%vis=VisibilityOfBellIneqWithNoChannel(PAPERbellcoeffs, LBoundPAPER_NUM, Pproj, CHSHstate, ins, outs)
p_entangled = ProbMultidimArray(CHSHstate, Pproj);
p_uniform = ProbMultidimArray(eye(4)/4, Pproj);
fprintf("b·p_ent=%f\n",sum(PAPERbellcoeffs .* p_entangled,'all'));
alpha = visibilityOfBellInequality(PAPERbellcoeffs, localbound, p_entangled, p_uniform);

%vis = 1;
%CHSHstate = vis * CHSHstate + (1-vis) * eye(size(CHSHstate))/size(CHSHstate,1);
ineqval = evaluate_bell_ineq(PAPERbellcoeffs, 0, eye(4)/4, Pproj)
1-alpha