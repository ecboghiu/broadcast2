basis = [1, "A0", "A1", "B0", "A0B0", "A1B0", "B1", "A0B1", "A1B1", ...
    "C0", "A0C0", "A1C0", "B0C0", "A0B0C0", "A1B0C0", "B1C0", "A0B1C0", "A1B1C0", ...
    "C1", "A0C1", "A1C1", "B0C1", "A0B0C1", "A1B0C1", "B1C1", "A0B1C1", "A1B1C1"];

basis_sym = [1];
for i=2:length(basis)
   %disp(basis(i))
    basis_sym = [basis_sym,  sym(basis(i))];
end

ineq = [-4 0 -2 0 0 0 -2 0 0 -1 -1 0 -1 0 -1 0 -1 1 -1 1 0 1 0 1 0 1 1];

A = [1,sym('A0'),sym('A1')];
B = [1,sym('B0'),sym('B1')];
C = [1,sym('C0'),sym('C1')];

basis_sym2 = kron(kron(C,B),A);


basis_sym2 - basis_sym

for 