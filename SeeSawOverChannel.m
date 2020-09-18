function [ChoiMap,optobj] = SeeSawOverChannel(alpha, bellcoeffs, Pproj, ins, outs)
    dimA = 2;
    dimB = 2;
    dimB1 = 2;
    dimB2 = 2;
    idB = eye(dimB);
    idB1 = eye(dimB1);
    idB2 = eye(dimB2);
    idB1B2 = Tensor(idB1,idB2);

    inputdimspace = 2;
    outputdimspace = 4;
    choidim = inputdimspace*outputdimspace;
    
    choi = sdpvar(choidim,choidim,'hermitian','complex');

    summ = 0;
    for x=ins{1}
        for y=ins{2}
            for z=ins{3}
                for a=outs{1}
                    for b=outs{2}
                        for c=outs{3}
                            term = Tensor(Pproj{1}{x}{a}, ...
                                            Pproj{2}{y}{b},...
                                            Pproj{3}{z}{c});
                            summ = summ + term*bellcoeffs(x,y,z,a,b,c);
                        end
                    end
                end
            end
        end        
    end
    
    state_small = ini_state(alpha);
    state_small_reshaped = reshape(state_small, 2,2,2,2);
    
    state = 0;
    for i=1:dimA
        for j=1:dimB
            for k=1:dimA
                for l=1:dimB
                    scnd = PartialTrace( choi * Tensor( ketbra(j,l,dimB).', idB1B2),...
                                        1, [dimB,dimB1*dimB2] );
                    state = state + state_small_reshaped(i,j,k,l) * ...
                                    Tensor(ketbra(i,k,dimA),scnd);
                end
            end
        end
    end
    
    constraints = [choi >= 0];
    constraints = [constraints, PartialTrace(choi, 2, [dimA, dimB1*dimB2]) == idB];
    
    objective = real(trace(state*summ));
    
    optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
    'verbose',0,'dualize',0, ...
    'showprogress',0,...
    'savesolverinput',1,'savesolveroutput',1,'debug',1));
    
    ChoiMap = value(choi);
    optobj = value(objective);
end
