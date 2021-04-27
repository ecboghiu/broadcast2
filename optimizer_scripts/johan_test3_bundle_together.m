%% TODO Johans suggestions to simplify the parameters

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
addpath(newdir);
addpath(newdir2);

%%

load('bellcoeffs_arxiv1112_2626.mat'); % loads 'bellcoeffs_cell','local_upper_bounds','ins','outs'
load('table3_arXiv1112_2626.mat'); % loads table3arXiv11122626

rng(123);

ineq_nr = 99;
bellcoeffs_ref = bellcoeffs_cell{ineq_nr};
localbound = local_upper_bounds(ineq_nr);

assert(mod(length(size(bellcoeffs_ref)),2)==0,"There should be as many inputs as outputs.");
dims = size(bellcoeffs_ref);
nrparties = length(dims)/2;
ins = dims(1:nrparties);
outs = dims(nrparties+1:end);

%% tests
a_ref = rand(4,4);
a_ref = (a_ref+a_ref')/2;
[a_r, a_i] = give_r_i(a_ref);
a = build_a(a_r, a_i);
assert(norm(a-a_ref)<1e-6,"func to take real and imag and go back doesnt work");



%% ini cond
povms = givePprojRANDgeneral(ins);
channel = {giveChannelRAND(2,4)};
bellcoeffs = bellcoeffs_ref;
state = NoisyWernerState(0);

tic
optimizer_ch = optimizer_Channel(ins, outs, 2);
toc
%%

bell_operator = give_Bell_operator(bellcoeffs, povms, ins, outs);

ia_povms = give_ia_povms(povms, ins, outs);
ia_state = give_ia_state(state);
ia_choi = give_ia_state(ChoiMatrix(channel));
ia_belloperator = give_ia_state(bell_operator);

% %%
% fprintf("Timing opt over channel");
% tic
% for i=1:25
%     [newchannel,newobjval1,~] = SeeSawOverChannel(state, bellcoeffs, povms);
% end
% toc
% fprintf('%f\n', newobjval1);
% 
% tic
% for i=1:25
%     output = optimizer_ch([{ia_state}, {bell_operator}]);
% end
% toc
% fprintf('%f\n', output{2});

%%
%partyidx = 2;
%tic
optimizer_alice = optimizer_povm_Alice(ins, outs);
optimizer_bob = optimizer_povm_Bob(ins, outs);



%toc
%optimizer_alice2 = optimizer_povm_party(partyidx, ins, outs);

%%
optimizer_a = optimizer_povm_party(1, ins, outs);
optimizer_b = optimizer_povm_party(2, ins, outs);
optimizer_c = optimizer_povm_party(3, ins, outs);

partial_products_for_a = give_partial_products(povms, bellcoeffs, 1, ins, outs);
partial_products_for_b = give_partial_products(povms, bellcoeffs, 2, ins, outs);
partial_products_for_c = give_partial_products(povms, bellcoeffs, 3, ins, outs);
%%

bellop1 = zeros(8,8);
for x=1:ins(1)
    for a=1:outs(1)
        bellop1 = bellop1 + kron(povms{1}{x}{a}, partial_products_for_a{x,a});
    end
end
bellop2 = zeros(8,8);
for y=1:ins(2)
    for b=1:outs(2)
        bellop2 = bellop2 + kron(povms{2}{y}{b}, partial_products_for_b{y,b});
    end
end
swapop=kron(SwapOperator([outs(2),outs(1)]),eye(outs(3)));
bellop2 = swapop * bellop2 * swapop';

bellop3 = zeros(8,8);
for z=1:ins(3)
    for c=1:outs(3)
        bellop3 = bellop3 + kron(partial_products_for_c{z,c}, povms{3}{z}{c});
    end
end

bellop = zeros(8,8);
for x=1:ins(1)
    for y=1:ins(2)
        for z=1:ins(3)
            for a=1:outs(1)
                for b=1:outs(2)
                    for c=1:outs(3)
                        bellop = bellop + bellcoeffs(x,y,z,a,b,c) * kron(kron(povms{1}{x}{a}, povms{2}{y}{b}), povms{3}{z}{c});
                    end
                end
            end
        end
    end
end
fprintf("Done");        

%%
fprintf("Timing opt over Alice's povm\n");

output_state = final_state(state, channel);
ia_output_state = give_ia_state(output_state);

[newPproj,newobjval2,~] = SeeSawOverASingleParty(2, output_state, bellcoeffs, povms);
output = optimizer_a([{ia_output_state}, partial_products_for_a(:)']);
fprintf("%f %f\n", newobjval2, output{1});

% fprintf("For a\n");
% tic
% for i=1:10
% [newPproj,newobjval2,~] = SeeSawOverASingleParty(1, output_state, bellcoeffs, povms);
% end
% toc
% fprintf('%f\n', newobjval2); 
% tic 
% for i=1:10
%     %output = optimizer_alice([{ia_output_state}, partial_products(:)']);
%     output = optimizer_a([{ia_output_state}, partial_products_for_a(:)']);
% end
% toc
% fprintf('%f\n', output{1}); 
% 
% 
% fprintf("For b\n");
% tic
% for i=1:10
% [newPproj,newobjval2,~] = SeeSawOverASingleParty(2, output_state, bellcoeffs, povms);
% end
% toc
% fprintf('%f\n', newobjval2); 
% tic 
% for i=1:10
%     %output = optimizer_alice([{ia_output_state}, partial_products(:)']);
%     output = optimizer_b([{ia_output_state}, partial_products_for_b(:)']);
% end
% toc
% fprintf('%f\n', output{1}); 
% 
% fprintf("For c\n");
% tic
% for i=1:10
% [newPproj,newobjval2,~] = SeeSawOverASingleParty(3, output_state, bellcoeffs, povms);
% end
% toc
% fprintf('%f\n', newobjval2); 
% tic 
% for i=1:10
%     %output = optimizer_alice([{ia_output_state}, partial_products(:)']);
%     output = optimizer_c([{ia_output_state}, partial_products_for_c(:)']);
% end
% toc
% fprintf('%f\n', output{1}); 

%%
% % fprintf("Timing over LP\n");
% % 
% % opt = optimizer_NS2_LP(ins, outs);
% % optimizerLP = opt{1};
% % nrconstr = opt{2};
% % 
% % rng(36);
% % povms = givePprojRANDgeneral(ins);
% % channel = {giveChannelRAND(2,4)};
% % bellcoeffs = bellcoeffs_ref;
% % state = NoisyWernerState(0);
% % 
% % fprintf("Optimizer NS2 LP");
% % tic
% % for i=1:10
% % p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), povms);
% % p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), povms);
% % [alpha, bellcoeffs, LPstatus, dual_alpha] = BroadcastNS2_LP(p_entangled, p_uniform);
% % end
% % toc
% % 
% % tic
% % for i=1:10
% % p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), povms);
% % p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), povms);
% % [alpha2, ~, ~, duals] = optimizerLP({p_entangled, p_uniform});
% % aux_shape = num2cell([ins, outs]);
% % bellcoeffs_opt = reshape(-duals{1}(1:prod([aux_shape{:}])), [aux_shape{:}]);
% % end
% % toc
% % 
% % fprintf("alpha=%f %f\n", alpha, alpha2);

%%

%%

%% funcs
function opt = optimizer_Channel(ins, outs, inputdimspace)
    fprintf("Calculating YALMIP's 'optimizer' of the scenario. It can take a while.\n");
    fprintf("Set 'OPTIMIZER_FLAG' to false to avoid 'optimizer' for debugging purposes.\n");
    nrparties = length(ins);

    dimA = outs(1);
    dimB = inputdimspace;
    dimB1 = outs(2);
    dimB2 = outs(3);
    idB = eye(dimB);
    idB1 = eye(dimB1);
    idB2 = eye(dimB2);
    idB1B2 = Tensor(idB1,idB2);

    %inputdimspace = dimB;
    outputdimspace = dimB1*dimB2;
    choidim = inputdimspace*outputdimspace;
    
    state_dim = dimA*dimB;
    state = sdpvar(state_dim,state_dim,'hermitian','complex');
    ia_state = give_ia_state(state);
    
%     aux_cell = num2cell(dims);
%     bellcoeffs = sdpvar(aux_cell{:},'full');
%     
%     povms = cell(nrparties, max(ins), max(outs));
%     ia_povms = cell(nrparties, max(ins), max(outs)); 
%     for party=1:nrparties
%         for in=1:ins(party)
%             for out=1:outs(party)
%                 povms{party,in,out} = sdpvar(outs(party),outs(party),'hermitian','complex');
%                 ia_povms{party,in,out} = povms{party,in,out}(find(triu(ones(outs(party)))));
%             end
%         end
%     end
    choi = sdpvar(choidim,choidim,'hermitian','complex');
    ia_choi = give_ia_state(choi);

    fprintf("Applying map through Choi isomorphism. \n");  
    state_small = state; %ini_state(alpha);
    state_small_reshaped = reshape(state_small, dimA,dimB,dimA,dimB);
    
    %%
   % TODO change this to a cleaner form, check that it gives the same
% %     state_sum = 0;
% %     for i=1:dimA
% %         for j=1:dimB
% %             for k=1:dimA
% %                 for l=1:dimB
% %                     %scnd = PartialTrace( choi * Tensor( ketbra(j,l,dimB).', idB1B2),...
% %                     %                    1, [dimB,dimB1*dimB2] );
% %                     state_sum = state_sum + state_small_reshaped(i,j,k,l) * ...
% %                                     kron(ketbra(i,k,dimA), ApplyMap(ketbra(j,l,dimB), choi));
% %                 end
% %             end
% %         end
% %     end
% %     output_state = state_sum;
    %%
%     biggerstate = kron( state.', eye(dimA*dimB1*dimB2) );
%     Phi = auxPHI(dimA);
%     biggerchannel = kron(Phi*Phi',choi);
%     % after previous line tensor spaces are A x A x B x B1 x B2, we need
%     % to swap the 2nd and 3rd
%     swapop=Tensor(eye(dimA),SwapOperator(2),eye(dimB1*dimB2)); % This ends up being a big matrix!
%     biggerchannel = swapop * biggerchannel * swapop';
%     output_state = PartialTrace(biggerchannel*biggerstate, [1,2], [dimA,dimB,dimA,dimB1,dimB2]);
%     %state = AppleMap(inistate, biggerchannel);  % either this or the
%                                                 % previous line work well
    %%
    biggerstate = Tensor(state.', eye(dimB1*dimB2) );
    biggerchannel = kron(eye(dimA), choi);
    output_state = PartialTrace(biggerchannel*biggerstate, 2, [dimA,dimB,dimB1,dimB2]);

                                                
    %%
    fprintf("Adding constraints. \n");  
    constraints_choi = cell(2);
    constraints_choi{1} = choi >= 0;
    constraints_choi{2} = PartialTrace(choi, 2, [dimB, dimB1*dimB2]) == idB;
%     
%     constraints_povm = cell(nrparties, max(ins));
%     for partyidx=1:nrparties
%         for x=1:ins(partyidx)
%             summ = 0;
%             for a=1:outs(partyidx)
%                 summ = summ + povms{party,x,a};
%                 constraints_povm{partyidx,x} = [constraints_povm{partyidx,x},  povms{party,x,a} >= 0];
%             end
%             constraints_povm{partyidx,x} = [constraints_povm{partyidx,x}, summ == eye(outs(partyidx))];
%         end
%     end
    
    fprintf("Calculating Bell operator.\n");
%     tic
%     state_times_BellOperator = 0;
%     for x=1:ins(1)
%         for y=1:ins(2)
%             for z=1:ins(3)
%                 for a=1:outs(1)
%                     for b=1:outs(2)
%                         for c=1:outs(3)
%                             fprintf("Progress on calculating symbolically output_state*BellOperator: %d %d %d %d %d %d\n", x, y, z, a, b, c);
%                             term = Tensor(povms{1,x,a}, ...
%                                             povms{2,y,b},...
%                                             povms{3,z,c});
%                             state_times_BellOperator = state_times_BellOperator + output_state*bellcoeffs(x,y,z,a,b,c)*term;
%                         end
%                     end
%                 end
%             end
%         end        
%     end
%     toc
    BellOperator = sdpvar(dimA*dimB1*dimB2,dimA*dimB1*dimB2,'full','complex');
    
    
    fprintf("Calculating the objective. \n");
    objective = real(trace(BellOperator*output_state));
    

    
    fprintf("Calculating optimizer object\n");
    opt_choi = optimizer([constraints_choi{:}], -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, {BellOperator}], ...
                        {choi, objective});
% 
%     opt_povm_A = optimizer([constraints_povm{1,:,:}], -objective, ...
%                          sdpsettings('solver','mosek'), ...
%                         [{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(2,:,:),1,[]),reshape(ia_povms(3,:,:),1,[])], ...
%                         [{objective}, povms(:)']);
%     opt_povm_B = optimizer([constraints_povm{2,:,:}], -objective, ...
%                          sdpsettings('solver','mosek'), ...
%                         [{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(1,:,:),1,[]),reshape(ia_povms(3,:,:),1,[])], ...
%                         [{objective}, povms(:)']);
%     opt_povm_C = optimizer([constraints_povm{3,:,:}], -objective, ...
%                          sdpsettings('solver','mosek'), ...
%                         [{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(1,:,:),1,[]),reshape(ia_povms(2,:,:),1,[])], ...
%                         [{objective}, povms(:)']);

    opt = opt_choi;
                    
% %     sol = optimize(constraints, -objective, ...
% %                     sdpsettings('solver','mosek','verbose', 0, ...
% %                                 'dualize', 0, 'showprogress', 0,...
% %                                 'savesolverinput', 0,'savesolveroutput', 0, ...
% %                                 'debug', 0, 'warning',0));
% %     if sol.problem ~= 0
% %        disp(sol);
% %        error('Check what problem there is.'); 
% %     end
% %     newChoiMap = value(choi);
% %     finalObj = value(objective);
% %     problemStatus = sol.problem;
    
    %list = whos;for i = 1:length(list);if strcmp(list(i).class,'sdpvar');clear(list(i).name);end;end
    %yalmip("clear");
end

function opt = optimizer_povm_Alice(ins, outs)
    fprintf("Calculating YALMIP's 'optimizer' of the scenario. It can take a while.\n");
    fprintf("Set 'OPTIMIZER_FLAG' to false to avoid 'optimizer' for debugging purposes.\n");
    nrparties = length(ins);

    dimA = outs(1);
    dimB1 = outs(2);
    dimB2 = outs(3);
    
    final_state_dim = dimA*dimB1*dimB2;
    final_state = sdpvar(final_state_dim,final_state_dim,'hermitian','complex');
    ia_state = give_ia_state(final_state);
    
    povms_alice = cell(ins(1), outs(1));
    partial_products = cell(ins(1), outs(1));
    %ia_partial_products = cell(ins(1), outs(1));
    for x=1:ins(1)
        for a=1:outs(1)
            povms_alice{x,a} = sdpvar(outs(1),outs(1),'hermitian','complex');
            dimen = outs(2)*outs(3);
            partial_products{x,a} = sdpvar(dimen,dimen,'full','complex');
            %ia_partial_products{x,a} = partial_products{x,a}(find(triu(ones(dimen)))); 
        end
    end

    partyidx = 1;
    constraints_povm = [];
    for x=1:ins(partyidx)
        summ = 0;
        for a=1:outs(partyidx)
            summ = summ + povms_alice{x,a};
            constraints_povm = [constraints_povm,  povms_alice{x,a} >= 0];
        end
        constraints_povm = [constraints_povm, summ == eye(outs(partyidx))];
    end

    
    BellOperator = 0;
    for x=1:ins(1)
       for a=1:outs(1)    
            BellOperator = BellOperator + Tensor(povms_alice{x,a}, partial_products{x,a});
       end
    end
    
   
    objective = real(trace(final_state*BellOperator));
    
    opt_povm_A = optimizer(constraints_povm, -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, partial_products(:)'], ...
                        [{objective}, povms_alice(:)']);

    opt = opt_povm_A;
end

function opt = optimizer_povm_Bob(ins, outs)
    fprintf("Calculating YALMIP's 'optimizer' of the scenario. It can take a while.\n");
    fprintf("Set 'OPTIMIZER_FLAG' to false to avoid 'optimizer' for debugging purposes.\n");
    nrparties = length(ins);

    dimA = outs(1);
    dimB1 = outs(2);
    dimB2 = outs(3);
    
    final_state_dim = dimA*dimB1*dimB2;
    final_state = sdpvar(final_state_dim,final_state_dim,'hermitian','complex');
    ia_state = give_ia_state(final_state);
    
    povms_bob = cell(ins(2), outs(2));
    partial_products = cell(ins(2), outs(2));
                dimen = outs(1)*outs(3);
    for y=1:ins(2)
        for b=1:outs(2)
            povms_bob{y,b} = sdpvar(outs(2),outs(2),'hermitian','complex');
            partial_products{y,b} = sdpvar(dimen,dimen,'full','complex');
        end
    end

    constraints_povm = [];
    for y=1:ins(2)
        summ = 0;
        for b=1:outs(2)
            summ = summ + povms_bob{y,b};
            constraints_povm = [constraints_povm,  povms_bob{y,b} >= 0];
        end
        constraints_povm = [constraints_povm, summ == eye(outs(2))];
    end

    BellOperator = 0;
    for y=1:ins(2)
       for b=1:outs(2)    
            BellOperator = BellOperator + Tensor(povms_bob{y,b}, partial_products{y,b});
       end
    end
    swapop = kron(SwapOperator([outs(2),outs(1)]),eye(outs(3)));
    BellOperator = swapop * BellOperator * swapop';
   
    objective = real(trace(final_state*BellOperator));
    
    opt_povm_A = optimizer(constraints_povm, -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, partial_products(:)'], ...
                        [{objective}, povms_bob(:)']);

    opt = opt_povm_A;
end

% TODO FINISH THIS
function opt = optimizer_povm_party(partyidx, ins, outs)
    fprintf("Calculating YALMIP's 'optimizer' of the scenario. It can take a while.\n");
    fprintf("Set 'OPTIMIZER_FLAG' to false to avoid 'optimizer' for debugging purposes.\n");
    nrparties = length(ins);

    dimA = outs(1);
    dimB1 = outs(2);
    dimB2 = outs(3);
    
    final_state_dim = dimA*dimB1*dimB2;
    final_state = sdpvar(final_state_dim,final_state_dim,'hermitian','complex');
    ia_state = give_ia_state(final_state);
    
    allbutone = logical(ones(nrparties,1));
    allbutone(partyidx) = false;
    prod_outs_without_partyidx = prod(outs(allbutone));
    
    povms_party = cell(ins(partyidx), outs(partyidx));
    partial_products = cell(ins(partyidx), outs(partyidx));
    %ia_partial_products = cell(ins(partyidx), outs(partyidx));
    for x=1:ins(partyidx)
        for a=1:outs(partyidx)
            povms_party{x,a} = sdpvar(outs(partyidx),outs(partyidx),'hermitian','complex');
            partial_products{x,a} = sdpvar(prod_outs_without_partyidx,prod_outs_without_partyidx,'full','complex');  % the combination with the bell coefficients can be arbitrary
            %ia_partial_products{x,a} = partial_products{x,a}(find(triu(ones(prod_outs_without_partyidx)))); 
        end
    end

    constraints_povm = [];
    for x=1:ins(partyidx)
        summ = 0;
        for a=1:outs(partyidx)
            summ = summ + povms_party{x,a};
            constraints_povm = [constraints_povm,  povms_party{x,a} >= 0];
        end
        constraints_povm = [constraints_povm, summ == eye(outs(partyidx))];
    end

    if partyidx == 2
        swapop = Tensor(SwapOperator([outs(2) outs(1)]), eye(outs(3)));
    end
    
    BellOperator = 0;
    for x=1:ins(partyidx)
       for a=1:outs(partyidx)
           if partyidx == 1  % Alice
               BellOperator = BellOperator + kron(povms_party{x,a}, partial_products{x,a});
           elseif partyidx == 2  % Bob
               % The idea is that we input 'partial_products' which
               % will be products of Alice times Charlie but we want to
               % tensor this with Bob. I do Bob x Alice x Charlie and then
               % use the 'swapop' to permute to Alice x Bob x Charlie
               BellOperator = BellOperator + swapop * kron(povms_party{x,a}, partial_products{x,a}) * swapop';
           elseif partyidx == 3  % Charlie
               BellOperator = BellOperator + kron(partial_products{x,a}, povms_party{x,a});
           else
              error('not supported for more than 3 parties'); 
           end
       end
    end
    
    objective = real(trace(final_state*BellOperator));
    
    opt_povm_party = optimizer(constraints_povm, -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, partial_products(:)'], ...
                        [{objective}, povms_party(:)']);

    opt = opt_povm_party;
end

function opt = optimizer_NS2_LP(ins, outs)
    nrparties = length(ins);
    
    nr_inputs_per_party = ins;
    nr_outputs_per_party = outs;
    
    auxdims = num2cell([ins,outs]);
    prob1 = sdpvar(auxdims{:},'full','real');
    prob2 = sdpvar(auxdims{:},'full','real'); 


    %nr_det_points = nroutputsofA^nrinputsofA;
    nr_det_points = nr_outputs_per_party.^nr_inputs_per_party;
  
    alpha = sdpvar(1); % alpha will be the visibility

    %% For lam
    tempdims = [nr_det_points(3), nr_inputs_per_party(1:2), nr_outputs_per_party(1:2)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    lamarray = sdpvar(prod(tempdims(:)),1);
    q_lam = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q_lam{coords{:}} = lamarray(idx);
    end
    
    %% For mu
    tempdims = [nr_det_points(2), nr_inputs_per_party(1), nr_inputs_per_party(3), nr_outputs_per_party(1), nr_outputs_per_party(3)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    muarray = sdpvar(prod(tempdims(:)),1);
    q_mu = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q_mu{coords{:}} = muarray(idx);
    end
    
    %% For nu
    tempdims = [nr_det_points(1), nr_inputs_per_party(2:end), nr_outputs_per_party(2:end)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    nuarray = sdpvar(prod(tempdims(:)),1);
    q_nu = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q_nu{coords{:}} = nuarray(idx);
    end
    

    visibility_constraints = [alpha >= 0];

    positivityconstraints = [];
	auxsize=size(q_lam);
    for i=1:prod(auxsize(:))
        positivityconstraints = [positivityconstraints, lamarray(i) >= 0];
    end
    auxsize=size(q_mu);
    for i=1:prod(auxsize(:))
        positivityconstraints = [positivityconstraints, muarray(i) >= 0];
    end
    auxsize=size(q_nu);
    for i=1:prod(auxsize(:))
        positivityconstraints = [positivityconstraints, nuarray(i) >= 0];
    end
    
    
    %% Non signalling constraints
    %% For q_nu
    nonsignalling_constraintsB_nu = [];
    % non signaling for bob:
    for nu = 1:nr_det_points(1)
        coordstructure = [nr_outputs_per_party(2), nr_inputs_per_party(2)];
        all_b_and_y = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_and_y,1)
            b = all_b_and_y(slice,1);
            y = all_b_and_y(slice,2);
            % choose z = 1
            summ1 = 0;
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q_nu{nu,y,1,b,c};
            end
            % the marginal for z != 1 should be equal to z=1
            for z=2:nr_inputs_per_party(3)
                summ2 = 0;
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q_nu{nu,y,z,b,c};
                end
                nonsignalling_constraintsB_nu = [nonsignalling_constraintsB_nu, summ1 == summ2];
            end
        end
    end
    %non signaling for charlie:
    nonsignalling_constraintsC_nu = [];
    for nu = 1:nr_det_points(1)
        coordstructure = [nr_outputs_per_party(3), nr_inputs_per_party(3)];
        all_b_and_y = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_and_y,1)
            c = all_b_and_y(slice,1);
            z = all_b_and_y(slice,2);
            % choose y = 1
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
                summ1 = summ1 + q_nu{nu,1,z,b,c};
            end
            % the marginal for z' != 1 should be equal to z=1
            for y=2:nr_inputs_per_party(2)
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                    summ2 = summ2 + q_nu{nu,y,z,b,c};
                end
                nonsignalling_constraintsC_nu = [nonsignalling_constraintsC_nu, summ1 == summ2];
            end
        end
    end
    % partial normalization nonsignaling (summing over bc should not depend
    % on the inputs)   
    nonsignalling_constraintsBC = [];
    for nu = 1:nr_det_points(1)          
        inputstructure = [nr_inputs_per_party(2), nr_inputs_per_party(3)];
        all_y_and_z = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        y1 = all_y_and_z(slice,1);
        z1 = all_y_and_z(slice,2);
        summ1 = 0;
        for b = 1:nr_outputs_per_party(2)
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q_nu{nu,y1,z1,b,c};
            end
        end
        for slice = 2:size(all_y_and_z,1)
            y2 = all_y_and_z(slice,1);
            z2 = all_y_and_z(slice,2);
            summ2 = 0;
            for b = 1:nr_outputs_per_party(2)
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q_nu{nu,y2,z2,b,c};
                end
            end
            nonsignalling_constraintsBC = [nonsignalling_constraintsBC, summ1 == summ2];
        end
    end
    
    %% For mu
    nonsignalling_constraintsA_mu = [];
    % non signaling for alice:
    for mu = 1:nr_det_points(2)
        coordstructure = [nr_outputs_per_party(1), nr_inputs_per_party(1)];
        all_a_and_x = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_a_and_x,1)
            a = all_a_and_x(slice,1);
            x = all_a_and_x(slice,2);
            % choose z = 1
            summ1 = 0;
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q_mu{mu,x,1,a,c};
            end
            % the marginal for z != 1 should be equal to z=1
            for z=2:nr_inputs_per_party(3)
                summ2 = 0;
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q_mu{mu,x,z,a,c};
                end
                nonsignalling_constraintsA_mu = [nonsignalling_constraintsA_mu, summ1 == summ2];
            end
        end
    end
    %non signaling for charlie:
    nonsignalling_constraintsC_mu = [];
    for mu = 1:nr_det_points(2)
        coordstructure = [nr_outputs_per_party(3), nr_inputs_per_party(3)];
        all_c_and_z = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_c_and_z,1)
            c = all_c_and_z(slice,1);
            z = all_c_and_z(slice,2);
            % choose x = 1
            summ1 = 0;
            for a = 1:nr_outputs_per_party(1)
                summ1 = summ1 + q_mu{mu,1,z,a,c};
            end
            % the marginal for z' != 1 should be equal to z=1
            for x=2:nr_inputs_per_party(1)
                summ2 = 0;
                for a = 1:nr_outputs_per_party(1)
                    summ2 = summ2 + q_mu{mu,x,z,a,c};
                end
                nonsignalling_constraintsC_mu = [nonsignalling_constraintsC_mu, summ1 == summ2];
            end
        end
    end
    % partial normalization nonsignaling (summing over bc should not depend
    % on the inputs)   
    nonsignalling_constraintsAC = [];
    for mu = 1:nr_det_points(2)          
        inputstructure = [nr_inputs_per_party(1), nr_inputs_per_party(3)];
        all_x_and_z = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        x1 = all_x_and_z(slice,1);
        z1 = all_x_and_z(slice,2);
        summ1 = 0;
        for a = 1:nr_outputs_per_party(1)
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q_mu{mu,x1,z1,a,c};
            end
        end
        
        for slice = 2:size(all_x_and_z,1)
            x2 = all_x_and_z(slice,1);
            z2 = all_x_and_z(slice,2);
            summ2 = 0;
            for a = 1:nr_outputs_per_party(1)
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q_mu{mu,x2,z2,a,c};
                end
            end
            nonsignalling_constraintsAC = [nonsignalling_constraintsAC, summ1 == summ2];
        end
    end
    
    %% For lam
    nonsignalling_constraintsA_lam = [];
    % non signaling for alice:
    for lam = 1:nr_det_points(3)
        coordstructure = [nr_outputs_per_party(1), nr_inputs_per_party(1)];
        all_a_and_x = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_a_and_x,1)
            a = all_a_and_x(slice,1);
            x = all_a_and_x(slice,2);
            % choose y = 1
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
                summ1 = summ1 + q_lam{lam,x,1,a,b};
            end
            % the marginal for y != 1 should be equal to y=1
            for y=2:nr_inputs_per_party(2)
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                    summ2 = summ2 + q_lam{lam,x,y,a,b};
                end
                nonsignalling_constraintsA_lam = [nonsignalling_constraintsA_lam, summ1 == summ2];
            end
        end
    end
    %non signaling for bob:
    nonsignalling_constraintsB_lam = [];
    for lam = 1:nr_det_points(3)
        coordstructure = [nr_outputs_per_party(2), nr_inputs_per_party(2)];
        all_b_and_y = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_and_y,1)
            b = all_b_and_y(slice,1);
            y = all_b_and_y(slice,2);
            % choose y = 1
            summ1 = 0;
            for a = 1:nr_outputs_per_party(1)
                summ1 = summ1 + q_lam{lam,1,y,a,b};
            end
            % the marginal for z' != 1 should be equal to z=1
            for x=2:nr_inputs_per_party(1)
                summ2 = 0;
                for a = 1:nr_outputs_per_party(1)
                    summ2 = summ2 + q_lam{lam,x,y,a,b};
                end
                nonsignalling_constraintsB_lam = [nonsignalling_constraintsB_lam, summ1 == summ2];
            end
        end
    end
    % partial normalization nonsignaling (summing over bc should not depend
    % on the inputs)   
    nonsignalling_constraintsAB = [];
    for lam = 1:nr_det_points(3)          
        inputstructure = [nr_inputs_per_party(1), nr_inputs_per_party(2)];
        all_x_and_y = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        x1 = all_x_and_y(slice,1);
        z1 = all_x_and_y(slice,2);
        summ1 = 0;
        for a = 1:nr_outputs_per_party(1)
            for b = 1:nr_outputs_per_party(2)
                summ1 = summ1 + q_lam{lam,x1,z1,a,b};
            end
        end
        
        for slice = 2:size(all_x_and_y,1)
            x2 = all_x_and_y(slice,1);
            z2 = all_x_and_y(slice,2);
            summ2 = 0;
            for a = 1:nr_outputs_per_party(1)
                for b = 1:nr_outputs_per_party(2)
                    summ2 = summ2 + q_lam{lam,x2,z2,a,b};
                end
            end
            nonsignalling_constraintsAB = [nonsignalling_constraintsAB, summ1 == summ2];
        end
    end
    
    
    %% Probability constraints
    det_strategyA = givedetstratA(nr_outputs_per_party(1),nr_inputs_per_party(1));
    det_strategyB = givedetstratA(nr_outputs_per_party(2),nr_inputs_per_party(2));
    det_strategyC = givedetstratA(nr_outputs_per_party(3),nr_inputs_per_party(3));

    coordstructure = [ins, outs];
    auxdims = num2cell([ins,outs]);
    noisy_prob = cell(auxdims{:});
    all_a_and_x = ind2subv(coordstructure, 1:prod(coordstructure(:)));
    for slice = 1:size(all_a_and_x,1)
        auxdims = num2cell(all_a_and_x(slice, :));
        noisy_prob{auxdims{:}} = (1-alpha) * prob1(auxdims{:}) + (alpha) * prob2(auxdims{:});
    end
    %noisy_prob = (1-alpha) * prob1 + (alpha) * prob2;
    
    probability_constraints = [];
    productstructure = [nr_inputs_per_party, nr_outputs_per_party];
    cartesianproduct_forprobconstraints = ind2subv(productstructure, 1:prod(productstructure(:)));
    for idx = 1:size(cartesianproduct_forprobconstraints,1)
        coords_cell = num2cell(cartesianproduct_forprobconstraints(idx,:));
        summ = 0;
        for nu = 1:nr_det_points(1)
            % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
            summ = summ + det_strategyA(nu, coords_cell{1}, coords_cell{nrparties+1}) ... 
                * q_nu{nu, coords_cell{2:nrparties}, coords_cell{nrparties+2:end}};
        end
        for mu = 1:nr_det_points(2)
            % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
            summ = summ + det_strategyB(mu, coords_cell{2}, coords_cell{nrparties+2}) ... 
                * q_mu{mu, coords_cell{1}, coords_cell{3}, coords_cell{nrparties+1}, coords_cell{nrparties+3}};
        end
        for lam = 1:nr_det_points(3)
            % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
            summ = summ + det_strategyC(lam, coords_cell{3}, coords_cell{nrparties+3}) ... 
                * q_lam{lam, coords_cell{1:2}, coords_cell{nrparties+1:nrparties+2}};
        end
        
        probability_constraints = [probability_constraints, summ == noisy_prob{coords_cell{:}}];
    end
    
    %% Solving the SDP
    
    objective = alpha;
            
    constraints = [probability_constraints, ...
        positivityconstraints, ...
        visibility_constraints, ...
        nonsignalling_constraintsA_lam, ...
        nonsignalling_constraintsA_mu, ...
        nonsignalling_constraintsB_lam, ...
        nonsignalling_constraintsB_nu, ...
        nonsignalling_constraintsC_nu, ...
        nonsignalling_constraintsC_mu, ...
        nonsignalling_constraintsBC, ...
        nonsignalling_constraintsAC, ...
        nonsignalling_constraintsAB];
    
    opt = cell(2);
    nr_duals = length(probability_constraints);
    opt{1} = optimizer(constraints, objective, ...
        sdpsettings('solver','mosek'), {prob1, prob2}, objective);
    opt{2} = nr_duals; % for extracting the LP
    
end

function [partial_products] = give_partial_products(povms, bellcoeffs, party_to_ignore, ins, outs)

nrparties = length(ins);
partial_products=cell(3);

if party_to_ignore == 1
    partial_products = cell(ins(1), outs(1));
 %   ia_partial_products = cell(ins(1), outs(1));
    dimen = outs(2)*outs(3);
    for x=1:ins(1)
        for a=1:outs(1)
        partial_products{x,a} = zeros(dimen,dimen);
            for y=1:ins(2)
                for z=1:ins(3)
                    for b=1:outs(2)
                        for c=1:outs(3)
                            partial_products{x,a} = partial_products{x,a} + bellcoeffs(x,y,z,a,b,c)*kron(povms{2}{y}{b},povms{3}{z}{c});                
                        end
                    end
                end
            end
 %       ia_partial_products{x,a} = partial_products{x,a}(find(triu(ones(dimen))));
        end
    end
elseif party_to_ignore == 2
% % %     partial_products = cell(ins(2), outs(2));
% % %  %   ia_partial_products = cell(ins(2), outs(2));
% % %     dimen = outs(1)*outs(3);
% % %     for y=1:ins(2)
% % %         for b=1:outs(2)
% % %         partial_products{y,b} = 0;
% % %             for x=1:ins(1)
% % %                 for z=1:ins(3)
% % %                     for a=1:outs(1)
% % %                         for c=1:outs(3)
% % %                             partial_products{y,b} = partial_products{y,b} + bellcoeffs(x,y,z,a,b,c)*kron(povms{1}{a}{x},povms{3}{z}{c});                
% % %                         end
% % %                     end
% % %                 end
% % %             end
% % %  %       ia_partial_products{y,b} = partial_products{y,b}(find(triu(ones(dimen))));
% % %         end
% % %     end
ppfb2 = cell(2,2);
for y=1:2
for b=1:2
ppfb2{y,b}=0;
for x=1:2
for a=1:2
for z=1:2
for c=1:2
ppfb2{y,b} = ppfb2{y,b} + bellcoeffs(x,y,z,a,b,c)*kron(povms{1}{x}{a},povms{3}{z}{c});
end
end
end
end
end
end
partial_products{2}=ppfb2;

elseif party_to_ignore == 3
    partial_products = cell(ins(3), outs(3));
 %   ia_partial_products = cell(ins(3), outs(3));
    dimen = outs(1)*outs(2);
    for z=1:ins(3)
        for c=1:outs(3)
        partial_products{z,c} = 0;
            for x=1:ins(1)
                for y=1:ins(2)
                    for a=1:outs(1)
                        for b=1:outs(2)
                            partial_products{z,c} = partial_products{z,c} + bellcoeffs(x,y,z,a,b,c)*kron(povms{1}{a}{x},povms{2}{y}{b});                
                        end
                    end
                end
            end
 %       ia_partial_products{z,c} = partial_products{z,c}(find(triu(ones(dimen))));
        end
    end

end
end

function [belloperator] = give_Bell_operator(bellcoeffs, povms, ins, outs)
belloperator = zeros(prod(outs(:)),prod(outs(:)));
for x=1:ins(1)
    for y=1:ins(2)
        for z=1:ins(3)
            for a=1:outs(1)
                for b=1:outs(2)
                    for c=1:outs(3)
                        term = Tensor(povms{1}{x}{a}, ...
                                        povms{2}{y}{b},...
                                        povms{3}{z}{c});
                        belloperator = belloperator + bellcoeffs(x,y,z,a,b,c)*term;
                    end
                end
            end
        end
    end        
end
end

function [a_r, a_i] = give_r_i(a)
n = size(a,1); 

indR =find(triu(ones(n)));
indI =find(triu(ones(n),1));

a_r =real(a(indR));
a_i =imag(a(indI));
end

function [a] = build_a(a_r, a_i)

n = (1 + sqrt(1+8*length(a_i)))/2;
assert(n==uint8(n), "not an integer, incorrect input!");
n = uint8(n);

indR =find(triu(ones(n)));
indI =find(triu(ones(n),1));

r1 = zeros(n,n);
r2 = zeros(n,n);
r1(indR) = a_r;
r2(indI) = a_i;
r1_nodiag = r1 - diag(diag(r1));
r1 = r1 + r1_nodiag';
r2 = r2 - r2';
a = r1 + 1j * r2;
end

function [ia_povms] = give_ia_povms(POVMS, ins, outs)

nrparties = length(ins);

ia_povms = cell(nrparties, max(ins), max(outs)); 
for party=1:nrparties
    for in=1:ins(party)
        for out=1:outs(party)
            ia_povms{party,in,out} = POVMS{party}{in}{out}(find(triu(ones(outs(party)))));
        end
    end
end
end

function [ia_state] = give_ia_state(state)
    ia_state = state(find(triu(ones(size(state,1)))));
end