function out = evaluate_bell_ineq(bellcoeffs, belloffset, finalstate, povms)
    %assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
    probarray = ProbMultidimArray(finalstate, povms);
	aux = probarray.*bellcoeffs;
    out = sum(aux(:)) + belloffset;
end

