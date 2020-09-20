function [newPproj,alpha] = SeeSawOverASingleParty(partyidx, state, bellcoeffs, Pproj, ins, outs)
    chosenPartyMeasurements = {{}};
    positivityconstraints  = [];
    dimofProjector = outs(partyidx);
    for x=1:ins(partyidx)
        for a=1:outs(partyidx)
            chosenPartyMeasurements{x}{a} = sdpvar(dimofProjector,dimofProjector,'hermitian','complex');
            positivityconstraints = [positivityconstraints, ...
                chosenPartyMeasurements{x}{a} >= 0];
            %sdisplay(chosenPartyMeasurements{x}{a} >= 0);
        end
    end
    
    % POVM constraints
    povmconstraints = [];
    for x=1:ins(partyidx)
        summ = 0;
        for a=1:outs(partyidx)
            summ = summ + chosenPartyMeasurements{x}{a};
        end
        povmconstraints = [povmconstraints, summ == eye(dimofProjector)];
        %sdisplay(summ == eye(dimofProjector));
    end

    summ = 0;
    for x=1:ins(1)
        for y=1:ins(2)
            for z=1:ins(3)
                for a=1:outs(1)
                    for b=1:outs(2)
                        for c=1:outs(3)
                            if partyidx == 1
                               term = Tensor(chosenPartyMeasurements{x}{a}, ...
                                            Pproj{2}{y}{b},...
                                            Pproj{3}{z}{c});
                            elseif partyidx == 2
                               term = Tensor(Pproj{1}{x}{a}, ...
                                            chosenPartyMeasurements{y}{b},...
                                            Pproj{3}{z}{c});
                            elseif partyidx == 3
                                term = Tensor(Pproj{1}{x}{a}, ...
                                            Pproj{2}{y}{b},...
                                            chosenPartyMeasurements{z}{c});
                            else
                                disp('error');
                            end
                            summ = summ + term*bellcoeffs(x,y,z,a,b,c);
                        end
                    end
                end
            end
        end        
    end
    
    objective = real(trace(state*summ));
    
    constraints = [positivityconstraints, povmconstraints];
    
    sol = optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
    'verbose',0,'dualize',0, ...
    'showprogress',0,...
    'savesolverinput',0,'savesolveroutput',0,'debug',0));

    newPproj = Pproj;
    for x=1:ins(partyidx)
        for a=1:outs(partyidx)
            newPproj{partyidx}{x}{a} = value(chosenPartyMeasurements{x}{a});
        end
    end
    %for x=ins{partyidx}
    %    fprintf("result for p=%d x=%d\n", partyidx, x);
    %    disp([newPproj{partyidx}{x}{1},newPproj{partyidx}{x}{2},newPproj{partyidx}{x}{1}+newPproj{partyidx}{x}{2}]); 
    %end
    
    alpha = value(objective);
end