rng('default');

ins = [4,4,2];
outs = [2,2,2];

C1 = sym('C1');
C2 = sym('C2');
B1 = sym('B1');
B2 = sym('B2');
B3 = sym('B3');
B4 = sym('B4');
A1 = sym('A1');
A2 = sym('A2');
A3 = sym('A3');
A4 = sym('A4');

EBI = A1*B1 + A1*B2 - A1*B3 - A1*B4 + A2*B1 - A2*B2 + ...
      A2*B3 - A2*B4 + A3*B1 - A3*B2 - A3*B3 + A3*B4;
LBoundEBI = 6;
OurEBI = EBI*(C1 + C2) + LBoundEBI*A4*(C1 - C2);
OurEBI = expand(OurEBI);

OurEBIInProb = ToProbabilityNotationIneqSym(OurEBI,ins,outs);
EBIbellcoeffs = GetBellCoeffsFromProbSymIneq(OurEBIInProb,ins,outs);

localMaximum = ClassicalOptInequality2(EBIbellcoeffs,ins,outs);
disp(localMaximum{1});