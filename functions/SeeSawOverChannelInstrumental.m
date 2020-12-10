function [newchannel,newobjval,problemStatus] = SeeSawOverChannelInstrumental(inistate, bellcoeffs, Pproj, ins_in, outs_in)
    dimA = 2;
    dimB = 2;
    dimB1 = 2;
    dimB2 = 2;
    idB = eye(dimB);
    idB1 = eye(dimB1);
    idB2 = eye(dimB2);
    idB1B2 = Tensor(idB1,idB2);

    %ins = size(bellcoeffs,1:3);
    %outs = size(bellcoeffs,5:7);
    %instr_ins = [size(bellcoeffs,4)];
    %instr_outs = [size(bellcoeffs,8)];
    ins = ins_in(1:3);
    outs = outs_in(1:3);
    instr_ins = ins_in(4);
    instr_outs = outs_in(4);
    
    
    inputdimspace = 2;
    outputdimspace = 4;
    choidim = inputdimspace*outputdimspace;
    
    choi = cell(instr_ins(1),instr_outs(1));
    constraints = [];
    for w = 1:instr_ins(1)
        for d = 1:instr_outs(1)
            choi{w}{d} = sdpvar(choidim,choidim,'hermitian','complex');
            constraints = [constraints, choi{w}{d} >= 0];
        end
    end
    
    for w = 1:instr_ins(1)
        summ = 0;
        for d = 1:instr_outs(1)
            summ = summ + choi{w}{d};
        end
        % after joining all the kraus operators we get a valid channel
        % which is a Choi map and which satisfies the following identity
        constraints = [constraints, PartialTrace(summ, 2, [dimB, dimB1*dimB2]) == idB];
    end

%     state_small = ini_state(alpha);
%     state_small_reshaped = reshape(state_small, 2,2,2,2);
    
    summ = 0;
    for w=1:instr_ins(1)
        for d=1:instr_outs(1)
    %         state = 0;
    %         for i=1:dimA
    %         for k=1:dimA
    %             
    %         for j=1:dimB
    %         for l=1:dimB
    %            % rho = rho_ijkl |ij><kl| 
    %            % sum_ijkl rho_ijkl |i><k| x tr_B [ Choi^w_d * ( |j><l|.T x Id_B1B2 )  ]
    %             scnd = PartialTrace( choi{w}{d} * Tensor( ketbra(j,l,dimB).', idB1B2),...
    %                                 1, [dimB,dimB1*dimB2] );
    %             state = state + state_small_reshaped(i,j,k,l) * ...
    %                             Tensor(ketbra(i,k,dimA),scnd);
    %         end
    %         end
    %         end
    %         end
            
            state_wd = ApplyMapToBipartition(inistate, choi{w}{d}, "right");
            
            summ_inner = 0;
            for x=1:ins(1)
                for y=1:ins(2)
                    for z=1:ins(3)
                        for a=1:outs(1)
                            for b=1:outs(2)
                                for c=1:outs(3)
                                    measurementOperators = Tensor(Pproj{1}{x}{a}, ...
                                                                  Pproj{2}{y}{b},...
                                                                  Pproj{3}{z}{c});
                                    summ_inner = summ_inner + bellcoeffs(x,y,z,w,a,b,c,d)*measurementOperators;
                                end
                            end
                        end
                    end
                end
            end
            summ_inner = summ_inner * state_wd; % just so I avoid multiplying by state_wd inside the loop, it should be slightly faster
            summ = summ + summ_inner;
        end
    end
    objective = real(trace(summ));
    

    
    sol = optimize(constraints, -objective, ...
                    sdpsettings('solver','mosek','verbose',0));
    channels = {{{}}};
    instr = 1; % NUMBER OF INSTRUMENTS
    for instrument=1:instr
       for w=1:instr_ins(1)
          for d=1:instr_outs(1)
              channels{instr}{w}{d} = value(choi{w}{d});
          end
       end
    end
    
    newchannel = channels;
    newobjval = value(objective);
    problemStatus = sol.problem;
end
