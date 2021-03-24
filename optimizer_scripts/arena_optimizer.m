load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads table3arXiv11122626

ineq_nr = 99;
bellcoeffs_ref = bellcoeffs_cell{ineq_nr};
localbound = local_upper_bounds(ineq_nr);

assert(mod(length(size(bellcoeffs_ref)),2)==0,"There should be as many inputs as outputs.");
dims = size(bellcoeffs_ref);
nrparties = length(dims)/2;
ins = dims(1:nrparties);
outs = dims(nrparties+1:end);

POVMs = givePprojRANDgeneral(ins);
restructPOVMs = {};
aux_idx = 1;
for party = 1:nrparties
   for input = 1:ins(party)
      for output = 1:outs(party)
          restructPOVMs{aux_idx} = POVMs{party}{input}{output};
          aux_idx = aux_idx + 1;
      end
   end
end
channel = {giveChannelRAND(2,4)};

%%
partyidx = 1;
%%

% p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs);
% p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), POVMs);

auxsize=num2cell([ins, outs]);
%coords = ind2subv(size(bellcoeffs), 1:prod(auxsize(:)));
bellcoeffs = sdpvar(auxsize{:});

%state_re = sdpvar(prod(outs(:)),prod(outs(:)),'full','real');%,'complex');
%state_im = sdpvar(prod(outs(:)),prod(outs(:)),'full','real');%,'complex');
%state = state_re + sqrt(-1) * state_im;
%state = (state + state')/2;
state = sdpvar(prod(outs(:)),prod(outs(:)),'hermitian','complex');%,'complex')

chosenPartyMeasurements = cell(ins(partyidx), outs(partyidx));
positivityconstraints  = [];
dimofProjector = outs(partyidx);
%chosenPartyMeasurements_split = {};
for x=1:ins(partyidx)
    for a=1:outs(partyidx)
 %       chosenP_re = sdpvar(dimofProjector,dimofProjector,'full','real');
 %       chosenP_im = sdpvar(dimofProjector,dimofProjector,'full','real');
 %       chosenP = chosenP_re + sqrt(-1) * chosenP_im;
 %       chosenPartyMeasurements_split{x, a, 1} = chosenP_re;
 %       chosenPartyMeasurements_split{x, a, 2} = chosenP_im;
        chosenPartyMeasurements{x,a} = sdpvar(dimofProjector,dimofProjector,'hermitian','complex');
        positivityconstraints = [positivityconstraints, chosenPartyMeasurements{x,a} >= 0];
    end
end


%% Define the SDP variables.
povms = cell(nrparties, max(ins), max(outs));
%povms_split = {};
for party = 1:nrparties
   for input = 1:ins(party)
      for output = 1:outs(party)
%          chosenP_re = sdpvar(dimofProjector,dimofProjector,'full','real');
%          chosenP_im = sdpvar(dimofProjector,dimofProjector,'full','real');
%          chosenP2 = chosenP_re + sqrt(-1) * chosenP_im;
%          povms_split{party,input,output,1} = chosenP_re;
%          povms_split{party,input,output,2} = chosenP_im;
          povms{party,input,output} = sdpvar(dimofProjector,dimofProjector,'hermitian','complex');
      end
   end
end
%% POVM constraints
povmconstraints = [];
for x=1:ins(partyidx)
    summ = 0;
    for a=1:outs(partyidx)
        summ = summ + chosenPartyMeasurements{x,a};
    end
    povmconstraints = [povmconstraints, summ == eye(dimofProjector)];
end

summ = 0;
for x=1:ins(1)
    for y=1:ins(2)
        for z=1:ins(3)
            for a=1:outs(1)
                for b=1:outs(2)
                    for c=1:outs(3)
                        if partyidx == 1
                           term = Tensor(chosenPartyMeasurements{x,a}, ...
                                        povms{2,y,b},...
                                        povms{3,z,c});
                        elseif partyidx == 2
                           term = Tensor(povms{1,x,a}, ...
                                        chosenPartyMeasurements{y,b},...
                                        povms{3,z,c});
                        elseif partyidx == 3
                            term = Tensor(povms{1,x,a}, ...
                                        povms{2,y,b},...
                                        chosenPartyMeasurements{z,c});
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
objective = -trace(state*summ);  % to be minimized, with -1 to max


%% Give the optimizer
opts = sdpsettings('solver','mosek', 'verbose', 2);

constraints = [positivityconstraints, povmconstraints];

params_in = {state, bellcoeffs, povms{1}}; 
solutions_out = {objective, chosenPartyMeasurements{1}};

optimizer_trial = optimizer(constraints, objective, opts, params_in, solutions_out);

num_state = final_state(NoisyWernerState(1-0.75), channel);
num_params_in = {num_state, bellcoeffs_ref, restructPOVMs{1}};%refactoredPOVMs{:}};
sol = optimizer_trial(num_params_in{:});
disp(sol{1});


[newPovms,finalAlpha,problemStatus] = SeeSawOverASingleParty(partyidx, num_state, bellcoeffs_ref, POVMs)
disp(finalAlpha);

%sol = optimize(constraints, -objective, opts);

% 
% % Return output
% newPovms = povms;
% for x=1:ins(partyidx)
%     for a=1:outs(partyidx)
%         newPovms{partyidx,x,a} = value(chosenPartyMeasurements{x}{a});
%     end
% end
% finalAlpha = value(objective);
% problemStatus = sol.problem;
%     
% list = whos;for i = 1:length(list);if strcmp(list(i).class,'sdpvar');clear(list(i).name);end;end
