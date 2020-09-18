function probbellexpression = ToProbabilityNotationIneqSym(corrbellexpression, ins, outs)
%ToProbabilityNotationTerm Summary of this function goes here
%   Detailed explanation goes here

expression = simplify(expand(corrbellexpression));

[C,T] = coeffs(expression);

nrterms = length(T);

summ = 0;
for termidx=1:nrterms
    summ = summ + C(termidx)*ToProbabilityNotationTerm(T(termidx),ins, outs);
end

probbellexpression = simplify(expand(summ));

end

