load bellcoeffs1.mat % load "bellcoefffs"
load symboliccorr1.mat % load "expression"

ins = [3,2,2];
outs = [2,2,2];


alpha = 0.57735;

C0 = sym('C0');
C1 = sym('C1');
B0 = sym('B0');
B1 = sym('B1');
A0 = sym('A0');
A1 = sym('A1');
A2 = sym('A2');
correct = (C0+C1) * (A1*(B0+B1) + A2*(B0-B1)) + 2*A0*(C0-C1);

% my conversion from p(abc|xyz) to <AxByCz>
[expression, scaling] = dispBellCoeffsCorrelators(bellcoeffs,ins, outs);
%disp(join([string(vpa(alpha*scaling,3)),'=',string(expression)]));

disp('correct:');
disp(expand(correct));

disp('mine:');
disp(expression);

disp('difference:');
diff = expand(correct)-expression;
% clean up small numbers
[C,T] = coeffs(diff);
tolerance = 1e-6;
C(abs(C)<tolerance)=0;
diff = vpa(dot(C,T),3); % vpa is to simplify integer fractions
disp(diff);