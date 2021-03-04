ins  = [3,2,2];
outs = [2,2,2];

insold = {[1,2,3],[1,2],[1,2]};
outsold = {[1,2],[1,2],[1,2]};

alpha = 0.57735;

C1 = sym('C1');
C2 = sym('C2');
B1 = sym('B1');
B2 = sym('B2');
A1 = sym('A1');
A2 = sym('A2');
A3 = sym('A3');

correct = (C1+C2) * (A2*(B1+B2) + A3*(B1-B2)) + 2*A1*(C1-C2);
correct = expand(correct);
disp(correct)

ab = ToProbabilityNotationIneqSym(correct,ins,outs);
disp(ab)

bellcoeffs = GetBellCoeffsFromProbSymIneq(ab,ins,outs);
%load bellcoeffs1.mat % load "bellcoefffs"
% load symboliccorr1.mat % load "expression
pexpression = dispBellCoeffs(bellcoeffs,insold,outsold);

bellcoeffs2 = GetBellCoeffsFromProbSymIneq(pexpression,ins,outs);

[expr,scaling] = dispBellCoeffsCorrelators(bellcoeffs,ins,outs);
disp(vpa(simplify(expr),3));

% conclusion: my function appears to be invertible