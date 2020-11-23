function out = ToCorrelatorNotationIneqSym(bellexpression,ins,outs)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

expression = simplify(expand(bellexpression));

[C,T] = coeffs(expression);

nrterms = length(T);

summ = 0;
for termidx=1:nrterms
    summ = summ + C(termidx)*ToCorrelatorNotationWithSymProbInput(T(termidx));
end

out = simplify(expand(summ));

end

