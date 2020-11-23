function [finalPproj,finalobj,finalchannel] = SeeSawOverAllPartiesInstrumental(bellcoeffs, initialState, initialPovms, initialChannels, ins, outs)
 
    ABS_TOL   = 1e-6;
    MAX_ITER  = 500;
    abschange = 1e6;
    objval    = -1e6;

    povms   = initialPovms;
    channels = initialChannels;
    
    party_ins = ins(1:3);
    party_outs = outs(1:3);
    instr_ins = [ins(4)];
    instr_outs = [outs(4)];

    %assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
    %nrparties = length(size(bellcoeffs))/2;
    nrparties = length(party_ins+1); % +1 for the channel
    
    partyidx = floor(1+(nrparties+1)*rand(1)); % start randomly
    iteration = 1;
    while abschange>ABS_TOL && iteration <= MAX_ITER
        if partyidx < nrparties+1
            [povms,newobjval,problemStatus] = SeeSawOverASinglePartyInstrumental(partyidx, initialState, bellcoeffs, povms, channels, ins, outs);
            abschange = abs(newobjval - objval);
            
            fprintf("iter=%d, partyidx=%d, bellineqvalue=%f, abschangeobjval=%g\n", iteration, partyidx, newobjval, abschange);
            objval = newobjval;
        else
            % now we do the channel
            [channels,newobjval,problemStatus] = SeeSawOverChannelInstrumental(initialState, bellcoeffs, povms, ins, outs);
                        
            abschange = abs(newobjval-objval);
            fprintf("iter=%d, channelq=%c, bellineqvalue=%f, abschangeobjval=%g\n", iteration, "~", newobjval, abschange);
   
            objval = newobjval;
        end
        partyidx = mod(partyidx,nrparties+1)+1;
        iteration = iteration + 1;
    end
    finalPproj = povms;
    finalobj = objval;
    finalchannel = channels; % FIX INSTR NUMBER
end