function output = SeeSawOverChannelInstr(alpha, bellcoeffs, Pproj, ins, outs)
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
    

    constraints = [];
    for w = 1:ins(4)
        for d = 1:outs(4)
            choi{w}{d} = sdpvar(choidim,choidim,'hermitian','complex');
            constraints = [constraints, choi{w}{d} >= 0];
        end
    end
    
    for w = 1:ins(4)
        summ = 0;
        for d = 1:outs(4)
            summ = summ + choi{w}{d};
        end
        % after joining all the kraus operators we get a valid channel
        % which is a Choi map and which satisfies the following identity
        constraints = [constraints, PartialTrace(summ, 2, [dimA, dimB1*dimB2]) == idB];
    end

    state_small = ini_state(alpha);
    state_small_reshaped = reshape(state_small, 2,2,2,2);
    
    summ = 0;
    for w=1:ins(4)
    for d=1:outs(4)
        state = 0;
        for i=1:dimA
        for j=1:dimB
        for k=1:dimA
        for l=1:dimB
            scnd = PartialTrace( choi{w}{d} * Tensor( ketbra(j,l,dimB).', idB1B2),...
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
                                measurementOperators = Tensor(Pproj{1}{x}{a}, ...
                                                              Pproj{2}{y}{b},...
                                                              Pproj{3}{z}{c});

                                summ = summ + bellcoeffs(x,y,z,w,a,b,c,d)*state*measurementOperators;
                            end
                        end
                    end
                end
            end
        end
    end        
    end
    objective = real(trace(summ));
    disp(objective)
    

    
    optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
    'verbose',0,'dualize',0, ...
    'showprogress',0,...
    'savesolverinput',1,'savesolveroutput',1,'debug',1));
    
    output = cell(1);
    output{1} = value(choi);
    output{2} = value(objective);
end
