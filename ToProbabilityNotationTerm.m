function probsym = ToProbabilityNotationTerm(symcorrelator, ins, outs)

factorarray = factor(symcorrelator);
nrparties = length(factorarray);

% I apologize but the following code is quite bad and ungeneralizable

if nrparties == 1 % something like Ax
    auxoutput  = GivePartyIdxAndInputFromSymName(factorarray);
    partyidx   = auxoutput(1);
    partyinput = auxoutput(2);
    
    if partyidx == 1 % can't think of a more elegant way to code it that using if elses
        x = partyinput;
        y = randi(ins(2)); % y,z don't matter because of no signaling, so I just choose them at random
        z = randi(ins(3));
        % now we must sum over b,c for different outputs of a
        prob1 = 0;
        prob2 = 0;
        for b=1:outs(2)
            for c=1:outs(3)
                prob1 = prob1 + GiveProbSymVarFromInputsOutputs([x,y,z],[1,b,c]);
                prob2 = prob2 + GiveProbSymVarFromInputsOutputs([x,y,z],[2,b,c]);
            end
        end     
        probsym = prob1 - prob2;
        
    elseif partyidx == 2
        x = randi(ins(1)); % x,z don't matter because of no signaling, so I just choose them at random
        y = partyinput;
        z = randi(ins(3));
        % now we must sum over b,c for different outputs of a
        prob1 = 0;
        prob2 = 0;
        for a=1:outs(1)
            for c=1:outs(3)
                prob1 = prob1 + GiveProbSymVarFromInputsOutputs([x,y,z],[a,1,c]);
                prob2 = prob2 + GiveProbSymVarFromInputsOutputs([x,y,z],[a,2,c]);
            end
        end
        probsym = prob1 - prob2;
    elseif partyidx == 3
        x = randi(ins(1)); % x,y don't matter because of no signaling, so I just choose them at random
        y = randi(ins(2));
        z = partyinput;
        % now we must sum over a,b for different outputs of a
        prob1 = 0;
        prob2 = 0;
        for a=1:outs(1)
            for b=1:outs(2)
                prob1 = prob1 + GiveProbSymVarFromInputsOutputs([x,y,z],[a,b,1]);
                prob2 = prob2 + GiveProbSymVarFromInputsOutputs([x,y,z],[a,b,2]);
            end
        end
        probsym = prob1 - prob2;
    else
        error("wrong");
    end
    
elseif nrparties == 2 % something like Ax*By
    % party1
    auxoutput  = GivePartyIdxAndInputFromSymName(factorarray(1));
    partyidx1   = auxoutput(1);
    partyinput1 = auxoutput(2);
    
    % party2
    auxoutput  = GivePartyIdxAndInputFromSymName(factorarray(2));
    partyidx2   = auxoutput(1);
    partyinput2 = auxoutput(2);
    
    if partyidx1 == partyidx2
       if partyinput1 == partyinput2
          probsym = ToProbabilityNotationTerm(factorarray(1),ins,outs);
          return;
       end
    end
    
    % again, really inelegant, but ill just do if else
    
    % first identify the missing party
    totnrofparties = length(outs);
    missingparty = setdiff(1:totnrofparties,[partyidx1,partyidx2]);
    
    if missingparty == 1
        x = randi(ins(1)); % take a random input 
        y = partyinput1;
        z = partyinput2;
        probsym = 0;
        for a=1:outs(missingparty)
           probsym = probsym +      GiveProbSymVarFromInputsOutputs([x,y,z],[a,1,1]) + ...
                         (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[a,1,2]) + ...
                         (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[a,2,1]) + ...
                              GiveProbSymVarFromInputsOutputs([x,y,z],[a,2,2]);
        end
        
    elseif missingparty == 2
        x = partyinput1;
        y = randi(ins(2)); % take a random input 
        z = partyinput2;
        probsym = 0;
        for b=1:outs(missingparty)
           probsym = probsym +      GiveProbSymVarFromInputsOutputs([x,y,z],[1,b,1]) + ...
                         (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[1,b,2]) + ...
                         (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[2,b,1]) + ...
                              GiveProbSymVarFromInputsOutputs([x,y,z],[2,b,2]);
        end
        
    elseif missingparty == 3
        x = partyinput1;
        y = partyinput2;
        z = randi(ins(3)); % take a random input 
        probsym = 0;
        for c=1:outs(missingparty)
           probsym = probsym +      GiveProbSymVarFromInputsOutputs([x,y,z],[1,1,c]) + ...
                         (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[1,2,c]) + ...
                         (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[2,1,c]) + ...
                              GiveProbSymVarFromInputsOutputs([x,y,z],[2,2,c]);
        end
        
    else
        error("shouldn't be here this code doesn't work for more than 3 parties") 
    end
    
    
    
elseif nrparties == 3 % something like Ax*By*Cz
    
    % party1
    auxoutput  = GivePartyIdxAndInputFromSymName(factorarray(1));
    partyidx1   = auxoutput(1);
    partyinput1 = auxoutput(2);
    
    % party2
    auxoutput  = GivePartyIdxAndInputFromSymName(factorarray(2));
    partyidx2   = auxoutput(1);
    partyinput2 = auxoutput(2);
    
    % party3
    auxoutput  = GivePartyIdxAndInputFromSymName(factorarray(3));
    partyidx3   = auxoutput(1);
    partyinput3 = auxoutput(2);
    
    % by default the ordering is Ax*By*Cz because its alphabetical, thus:
    x=partyinput1;
    y=partyinput2;
    z=partyinput3;
    
    probsym =    GiveProbSymVarFromInputsOutputs([x,y,z],[1,1,1]) + ...
            (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[1,2,1]) + ...
            (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[2,1,1]) + ...
                 GiveProbSymVarFromInputsOutputs([x,y,z],[2,2,1]) + ...
            (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[1,1,2]) + ...
                 GiveProbSymVarFromInputsOutputs([x,y,z],[1,2,2]) + ...
                 GiveProbSymVarFromInputsOutputs([x,y,z],[2,1,2]) + ...
            (-1)*GiveProbSymVarFromInputsOutputs([x,y,z],[2,2,2]);
    else
    error("Shouldn't be at this point");
end
    

end

function out = GivePartyIdxAndInputFromSymName(sym)
% input something like A0 or B1 or C0
% example: input B1
% output: [2,1], substitute A->1, B->2, C->3
name = char(sym); % ex: sym var A0 --> 'A0'
partyIdx = name(1)-'A'+1;
partyInput = str2double(name(2));
out(1) = partyIdx;
out(2) = partyInput;
end

function symvar = GiveProbSymVarFromInputsOutputs(settings,outputs)
x=settings(1);
y=settings(2);
z=settings(3);
a=outputs(1);
b=outputs(2);
c=outputs(3);
str = join(['p_',string(a),...
                 string(b),...
                 string(c),...
                 '_',...
                 string(x),...
                 string(y),...
                 string(z)],'');
symvar = sym(char(str));
end