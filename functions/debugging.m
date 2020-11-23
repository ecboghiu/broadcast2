ins = {[1,2,3,4],[1,2],[1,2,3,4]};
outs = {[1,2],[1,2],[1,2]};

A_0 = sym('A_0');
A_1 = sym('A_1');
A_2 = sym('A_2');
A_3 = sym('A_3');
A_4 = sym('A_4');

B_0 = sym('B_0');
B_1 = sym('B_1');

C_0 = sym('C_0');
C_1 = sym('C_1');
C_2 = sym('C_2');
C_3 = sym('C_3');


EBI = A_0*C_0 + A_0*C_1 - A_0*C_2 - A_0*C_3 + A_1*C_0 - A_1*C_1 + ...
      A_1*C_2 - A_1*C_3 + A_2*C_0 - A_2*C_1 - A_2*C_2 + A_2*C_3;

I = EBI*(B_1 + B_0) + 6*A_5*(B_1 - B_0);
I = expand(I);

ins = {[1,2,3],[1,2],[1,2]};
outs = {[1,2],[1,2],[1,2]};
out = ClassicalOptInequality(bellcoeffs, ins, outs);
disp(out{1})
disp(out{2})

