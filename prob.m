function p = prob(state,Pproj,outputvec,inputvec)
    nrparties = length(Pproj);
    tensor = 1;
    for idx=1:nrparties
        inp = inputvec(idx);
        out = outputvec(idx);
        tensor = kron(tensor,Pproj{idx}{inp}{out});
    end

    p = real(trace(state*tensor));
end