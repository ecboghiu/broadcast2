function out = cleanChannel(channel, dim_in, dim_out, tol)
tol_psd = tol;
tol_cleansmallnumbers = 1e-12;

% input must be in choi form
Cm = clean( ChoiMatrix(channel) , tol_cleansmallnumbers);
% Cm = (Cm + Cm')/2; % make hermitian

d1 = dim_in;
d2 = dim_out;

two_C = kron(PartialTrace(Cm, 2, [d1, d2]), eye(d2)/d2);
onetwo_C = trace(Cm)*eye(d1*d2)/(d1*d2);

Cprime = Cm - two_C + onetwo_C;
Cprimeprime = Cprime*d1/trace(Cm);

[eigvec,eigvals] = eig(Cprimeprime);

smallesteig = min(clean(diag(eigvals),tol_cleansmallnumbers));
if smallesteig < 0
    psdCdiag = (1-abs(smallesteig)) * diag(eigvals) + abs(smallesteig) * eye(d1*d2)/d2;
    out = eigvec^-1 * psdCdiag * eigvec;
else
    out = Cprimeprime;
end

end

