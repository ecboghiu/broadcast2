B1 = sym('B1');
B2 = sym('B2');
B3 = sym('B3');
B4 = sym('B4');
A1 = sym('A1');
A2 = sym('A2');
A3 = sym('A3');

ins = [3,4];
outs = [2,2];

EBI = A1*B1 + A1*B2 - A1*B3 - A1*B4 + A2*B1 - A2*B2 + ...
      A2*B3 - A2*B4 + A3*B1 - A3*B2 - A3*B3 + A3*B4;
%   
% EBI2 = B1*A1 + B1*A2 - B1*A3 - B1*A4 + B2*A1 - B2*A2 + ...
%       B2*A3 - B2*A4 + B3*A1 - B3*A2 - B3*A3 + B3*A4;

EBInProb = ToProbabilityNotationIneqSym(EBI,ins,outs);
EBIbellcoeffs = GetBellCoeffsFromProbSymIneq(EBInProb,ins,outs);
%LBoundOurEBI = ClassicalOptInequality2(OurEBI,ins,outs);
%LBoundOurEBI = ClassicalOptInequality_fromLPBroadcast(EBIbellcoeffs);
LBoundOurEBI = ClassicalOptInequalityNoBroadcast(EBIbellcoeffs);
fprintf('LBoundOurEBI: %f\n' ,LBoundOurEBI);

maxentstate = 1/sqrt(2)*[1 0 0 1].';
maxentstate = maxentstate * maxentstate';

sig1 = [0+0i,1+0i;1+0i,0+0i];
sig2 = [0+0i,-1i;1i,0+0i];
sig3 = [1+0i,0+0i;0+0i,-1+0i];

nA0 = [0,0,1]; % get sigma_Z
nA1 = [1,0,0]; % get sigma_X
nA2 = [0,1,0]; % get sigma_Y
A0 = nA0(1) * sig1 + nA0(2) * sig2 + nA0(3) * sig3;
A1 = nA1(1) * sig1 + nA1(2) * sig2 + nA1(3) * sig3;
A2 = nA2(1) * sig1 + nA2(2) * sig2 + nA2(3) * sig3;

nB0 = 1/sqrt(3) * [1,-1,1];
nB1 = 1/sqrt(3) * [-1,1,1];
nB2 = 1/sqrt(3) * [1,1,-1];
nB3 = 1/sqrt(3) * [-1,-1,-1];
B0 = nB0(1) * sig1 + nB0(2) * sig2 + nB0(3) * sig3;
B1 = nB1(1) * sig1 + nB1(2) * sig2 + nB1(3) * sig3;
B2 = nB2(1) * sig1 + nB2(2) * sig2 + nB2(3) * sig3;
B3 = nB3(1) * sig1 + nB3(2) * sig2 + nB3(3) * sig3;

povms = {{giveprojs(A0),giveprojs(A1),giveprojs(A2)}, ...
            {giveprojs(B0),giveprojs(B1),giveprojs(B2),giveprojs(B3)}};
        
p_entangled = ProbMultidimArray(maxentstate, povms);
p_uniform = ProbMultidimArray(eye(4)/4, povms);

fprintf("This should be 6.9282");
ineqval=evaluate_bell_ineq(EBIbellcoeffs, 0, maxentstate, povms)
evaluate_bell_ineq(EBIbellcoeffs, 0, eye(4)/4, povms)

alpha = visibilityOfBellInequality(EBIbellcoeffs, LBoundOurEBI, p_entangled, p_uniform)



