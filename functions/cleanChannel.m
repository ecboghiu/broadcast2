function out = cleanChannel(channel, d1, d2)
if ~IsPSD(channel)
    Cm = ChoiMatrix(channel);  % input must be in choi form

    % d1 = in
    % d2 = out
    two_C = kron(PartialTrace(Cm, 2, [d1, d2]), eye(d2)/d2);
    onetwo_C = trace(Cm)*eye(d1*d2)/(d1*d2);

    Cprime = Cm - two_C + onetwo_C;
    Cprimeprime = Cprime*d1/trace(Cprime);

    smallesteig = min(eig(Cprimeprime));
    if smallesteig < 0
        eta = abs(smallesteig)/(abs(smallesteig)+1/d2);
        out = (1-eta) * Cprimeprime + eta * eye(d1*d2)/d2;
    else
        out = Cprimeprime;
    end

    assert(norm(PartialTrace(out,2,[d1,d2])-eye(d1),'fro')<1e-12,"Choi matrix condition not satisfied");
    assert(IsPSD(out),"Choi matrix not positive");
    fprintf("after cleaning: %f smallest_eig=%f\n", norm(out-channel, 'fro'), smallesteig);
    
else
    out=channel;
end

end

