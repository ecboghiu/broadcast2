ins = [3,4,2];
outs = [2,2,2];


A1 = sym('A1');
A2 = sym('A2');
A3 = sym('A3');
B1 = sym('B1');
B2 = sym('B2');
B3 = sym('B3');
B4 = sym('B4');

%So that the program works, but the C's won't appear anywhere
C1 = sym('C1');
C2 = sym('C2');

EBI = A1*B1 + A1*B2 - A1*B3 - A1*B4 + A2*B1 - A2*B2 + ...
      A2*B3 - A2*B4 + A3*B1 - A3*B2 - A3*B3 + A3*B4;
  
EBIinProb = ToProbabilityNotationIneqSym(EBI,ins,outs);
EBIbellcoeffs = GetBellCoeffsFromProbSymIneq(EBIinProb,ins,outs);
LBoundEBI = ClassicalOptInequality_fromLPBroadcast(EBIbellcoeffs);
fprintf('LBoundOurEBI: %f\n',LBoundEBI);