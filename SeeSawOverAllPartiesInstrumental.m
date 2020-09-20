function [finalPproj,finalobj,finalchannel] = SeeSawOverAllParties(inistate, bellcoeffs, ...
                                            initialPproj, ins, outs, alpha, inichannel)
    nrparties = length(ins);
    objval = -1e6;

    Pproj = initialPproj;
    ABS_TOL = 1e-8;
    MAX_ITER = 1000;
    iteration = 1;
    abschange = 1e6;
    partyidx = 1;
    channel = inichannel;
    state = final_state(ini_state(alpha),channel);%inistate;
    
    dims = num2cell([ins,outs]);
    probability = zeros(dims{:}); % call as probability(x,y,z,a,b,c) for 3 parties

    while abschange>ABS_TOL && iteration <= MAX_ITER
        if partyidx < nrparties+1
            [newPproj,newobjval] = SeeSawOverASingleParty(partyidx,state, bellcoeffs, Pproj, ins, outs);
            Pproj = newPproj;
            
            abschange = abs(newobjval - objval);
            
            %fprintf("iter=%d, partyidx=%d, bellineqvalue=%f, abschangeobjval=%g\n", iteration, partyidx, newobjval, abschange);
            
            objval = newobjval;
        else
            % now we do the channel
            
            [newchannel,newobjval] = SeeSawOverChannel(alpha, bellcoeffs, Pproj, ins, outs);
            channel = newchannel;
            state = final_state(ini_state(alpha),channel);
            
            
    for x=1:ins(1)
        for y=1:ins(2)
            for z=1:ins(3)
                for a=1:outs(1)
                    for b=1:outs(2)
                        for c=1:outs(3)
            probability(x,y,z,a,b,c)=prob(state,Pproj,[a,b,c],[x,y,z]);
                        end 
                    end 
                end 
            end 
        end 
    end 
            checkThatProbSumsToOne(probability, ins, outs);
            
            abschange = abs(newobjval-objval);
            %fprintf("iter=%d, channelq=%c, bellineqvalue=%f, abschangeobjval=%g\n", iteration, "~", newobjval, abschange);
            objval = newobjval;
        end
        partyidx = mod(partyidx,nrparties+1)+1;
        iteration = iteration + 1;
    end
    finalPproj = Pproj;
    finalobj = objval;
    finalchannel = channel;
end