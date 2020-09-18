function bellcoeffs = GetBellCoeffsFromProbSymIneq(symineq,ins,outs)
% WARNING this only works for 3 parties

dims = num2cell([ins,outs]);
bellcoeffs = zeros(dims{:}); % call as bellcoeffs(x,y,z,a,b,c) for 3 parties

[C,T] = coeffs(symineq);
nrterms = length(T);

for i=1:nrterms
    nameinput = char(T(i));
    % assuming input of the form p_abc_xyz and no integer greater than 10

    a = str2double(nameinput(3));
    b = str2double(nameinput(4));
    c = str2double(nameinput(5));
    x = str2double(nameinput(7));
    y = str2double(nameinput(8));
    z = str2double(nameinput(9));
       
    bellcoeffs(x,y,z,a,b,c)=C(i);
end


end

