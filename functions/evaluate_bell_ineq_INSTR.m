function out = evaluate_bell_ineq_INSTR(bellcoeffs, belloffset, finalstate, povms, channels)
    assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
    probarray = ProbMultidimArrayInstrumental(finalstate,povms,channels);
    out = sum(probarray.*bellcoeffs,'all') + belloffset;
end
