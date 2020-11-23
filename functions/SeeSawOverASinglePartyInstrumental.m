function [newPproj,bellineqvalue,problemStatus] = SeeSawOverASinglePartyInstrumental(partyidx, initialState, bellcoeffs, Pproj, channels, ins_in, outs_in)
    dimA = 2;
    dimB = 2;
    dimB1 = 2;
    dimB2 = 2;
    idB = eye(dimB);
    idB1 = eye(dimB1);
    idB2 = eye(dimB2);
    idB1B2 = Tensor(idB1,idB2);

    % BY HAND CHANGE LATER
    %ins = size(bellcoeffs,1:3);
    %outs = size(bellcoeffs,5:7);
    %instr_ins = [size(bellcoeffs,4)];
    %instr_outs = [size(bellcoeffs,8)];
    ins = ins_in(1:3);
    outs = outs_in(1:3);
    instr_ins = ins_in(4);
    instr_outs = outs_in(4);
    
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

    state_small = initialState;
    state_small_reshaped = reshape(state_small, 2,2,2,2);
    instr=1;
    summ = 0;
    for w=1:instr_ins(1)
    for d=1:instr_outs(1)
        state = 0;
        for i=1:dimA
        for j=1:dimB
        for k=1:dimA
        for l=1:dimB
            scnd = PartialTrace( ChoiMatrix(channels{instr}{w}{d}) * Tensor( ketbra(j,l,dimB).', idB1B2),...
                                1, [dimB,dimB1*dimB2] );
            state = state + state_small_reshaped(i,j,k,l) * ...
                            Tensor(ketbra(i,k,dimA),scnd);
        end
        end
        end
        end
        for x=1:ins(1)
            for y=1:ins(2)
                for z=1:ins(3)
                    for a=1:outs(1)
                        for b=1:outs(2)
                            for c=1:outs(3)
                                if partyidx == 1
                                   measterm = Tensor(chosenPartyMeasurements{x}{a}, ...
                                                Pproj{2}{y}{b},...
                                                Pproj{3}{z}{c});
                                elseif partyidx == 2
                                   measterm = Tensor(Pproj{1}{x}{a}, ...
                                                chosenPartyMeasurements{y}{b},...
                                                Pproj{3}{z}{c});
                                elseif partyidx == 3
                                    measterm = Tensor(Pproj{1}{x}{a}, ...
                                                Pproj{2}{y}{b},...
                                                chosenPartyMeasurements{z}{c});
                                else
                                    disp('error');
                                end
                                summ = summ + bellcoeffs(x,y,z,w,a,b,c,d)*measterm*state;
                            end
                        end
                    end
                end
            end        
        end
    end
    end
    
    objective = real(trace(summ));
    
    constraints = [positivityconstraints, povmconstraints];
    
    sol = optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
    'verbose',0,'dualize',0, ...
    'showprogress',0,...
    'savesolverinput',0,'savesolveroutput',0,'debug',0));
    problemStatus = sol.problem;
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
    
    bellineqvalue = value(objective);
end