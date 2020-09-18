function out = MetaCompleteRoundFixedBellCorrIneqSym(ineq, ins, outs)

OurEBIInProb = ToProbabilityNotationIneqSym(OurEBI,ins,outs);
EBIbellcoeffs = GetBellCoeffsFromProbSymIneq(OurEBIInProb,ins,outs);

LBoundOurEBI = ClassicalOptInequality2(EBIbellcoeffs,ins,outs);


% localMaximum = ClassicalOptInequality(EBIbellcoeffs,ins,outs);
% disp(localMaximum{1});
% disp(localMaximum{2});

output = CompleteRoundFixedBell(EBIbellcoeffs, LBoundOurEBI, ins, outs);
alpha = output{1};

TOL = 1e-1;
MAXITER = 1000;

iteration = 1;
while ( abs(alpha-1)<TOL || abs(alpha-0)<TOL || alpha-0.707>0 || alpha-0.577>0 ) && iteration <= MAXITER
    output =  CompleteRoundFixedBell(EBIbellcoeffs, LBoundOurEBI, ins, outs);
    alpha = output{1};
    fprintf("\n New meta iter. \n\n");
    
    iteration = iteration + 1;
end
    
out = output;
end


