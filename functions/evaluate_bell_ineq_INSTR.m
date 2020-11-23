function out = evaluate_bell_ineq_INSTR(bellcoeffs, belloffset, finalstate, povms, channels)
    %assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");  % remove this condition because we might have trailing 1-dims in the vector which gets annoying
    probarray = ProbMultidimArrayInstrumental(finalstate,povms,channels);
	aux=probarray.*bellcoeffs;
    out = sum(aux(:)) + belloffset;
end
