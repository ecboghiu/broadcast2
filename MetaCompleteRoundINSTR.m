function out = MetaCompleteRoundINSTR(ins, outs)

output = CompleteRoundInstr(0,0,ins,outs);
alpha = output{1};

TOL = 1e-1;
MAXITER = 100;

iteration = 1;
while ( abs(alpha-1)<TOL || abs(alpha-0)<TOL || alpha-0.707>0 || alpha-0.576>0) && iteration <= MAXITER
    output = CompleteRoundInstr(0,0,ins,outs);
    alpha = output{1};
    fprintf("\n New meta iter. \n\n");
    
    iteration = iteration + 1;
end
    
out = output;
end

