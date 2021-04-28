mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
addpath(newdir);
addpath(newdir2);

%%
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


bell_operator = give_Bell_operator(bellcoeffs, povms, ins, outs);
ia_povms = give_ia_povms(povms, ins, outs);
ia_state = give_ia_state(state);
ia_choi = give_ia_state(ChoiMatrix(channel));
ia_belloperator = give_ia_state(bell_operator);

%%

tic
optimizer_ch = optimizer_channel(ins, outs, 2);
toc

fprintf("Loop time improvement for using 'optimizer' with channel.\n");
tic
for i=1:25
    [newchannel,newobjval1,~] = SeeSawOverChannel(state, bellcoeffs, povms);
end
toc
fprintf('%f\n', newobjval1);

tic
for i=1:25
    output = optimizer_ch([{ia_state}, {bell_operator}]);
end
toc
fprintf('%f\n', output{2});

assert(abs(newobjval1-output{2}) <1e-10, "There should be the same output");

%%
optimizer_a = optimizer_povm_party(1, ins, outs);
optimizer_b = optimizer_povm_party(2, ins, outs);
optimizer_c = optimizer_povm_party(3, ins, outs);

partial_products_for_a = give_partial_products(povms, bellcoeffs, 1, ins, outs);
partial_products_for_b = give_partial_products(povms, bellcoeffs, 2, ins, outs);
partial_products_for_c = give_partial_products(povms, bellcoeffs, 3, ins, outs);   

%%
fprintf("Loop time improvement for using 'optimizer' with the POVM optimization.\n");

output_state = final_state(state, channel);
ia_output_state = give_ia_state(output_state);

fprintf("For Alice\n");
tic
for i=1:10
[newPproj,newobjval2,~] = SeeSawOverASingleParty(1, output_state, bellcoeffs, povms);
end
toc
fprintf('%f\n', newobjval2); 
tic 
for i=1:10
    output = optimizer_a([{ia_output_state}, partial_products_for_a(:)']);
    povm_alice = reshape({output{2:end}},[ins(1),outs(1)]);
    newPproj2 = newPproj;
    for x=1:ins(1)
        for a=1:outs(1)
            newPproj2{1}{x}{a} = povm_alice{x,a};
        end
    end
    for party=1:nrparties
        for x=1:ins(party)
            for a=1:outs(party)
                assert(norm(newPproj{party}{x}{a}-newPproj2{party}{x}{a})<1e-10,"Output povms are not the same with both methods");
            end
        end
    end
end
toc
fprintf('%f\n', output{1}); 
assert(abs(newobjval2-output{1}) <1e-12, "There should be the same output");

fprintf("For Bob\n");
tic
for i=1:10
[newPproj,newobjval2,~] = SeeSawOverASingleParty(2, output_state, bellcoeffs, povms);
end
toc
fprintf('%f\n', newobjval2); 
tic 
for i=1:10
    output = optimizer_b([{ia_output_state}, partial_products_for_b(:)']);
    povm_bob = reshape({output{2:end}},[ins(2),outs(2)]);
    newPproj2 = newPproj;
    for x=1:ins(2)
        for a=1:outs(2)
            newPproj2{2}{x}{a} = povm_bob{x,a};
        end
    end
    for party=1:nrparties
        for x=1:ins(party)
            for a=1:outs(party)
                assert(norm(newPproj{party}{x}{a}-newPproj2{party}{x}{a})<1e-10,"Output povms are not the same with both methods");
            end
        end
    end
end
toc
fprintf('%f\n', output{1}); 
assert(abs(newobjval2-output{1}) <1e-12, "There should be the same output");

fprintf("For Charlie\n");
tic
for i=1:10
[newPproj,newobjval2,~] = SeeSawOverASingleParty(3, output_state, bellcoeffs, povms);
end
toc
fprintf('%f\n', newobjval2); 
tic 
for i=1:10
    output = optimizer_c([{ia_output_state}, partial_products_for_c(:)']);
    povm_charlie = reshape({output{2:end}},[ins(3),outs(3)]);
    newPproj2 = newPproj;
    for x=1:ins(3)
        for a=1:outs(3)
            newPproj2{3}{x}{a} = povm_charlie{x,a};
        end
    end
    for party=1:nrparties
        for x=1:ins(party)
            for a=1:outs(party)
                assert(norm(newPproj{party}{x}{a}-newPproj2{party}{x}{a})<1e-10,"Output povms are not the same with both methods");
            end
        end
    end
end
toc
fprintf('%f\n', output{1}); 
assert(abs(newobjval2-output{1}) <1e-12, "There should be the same output");

%%
fprintf("Timing over LP\n");

opt = optimizer_NS2_LP(ins, outs);
optimizerLP = opt{1};
nrconstr = opt{2};

rng(36);
povms = givePprojRANDgeneral(ins);
channel = {giveChannelRAND(2,4)};
bellcoeffs = bellcoeffs_ref;
state = NoisyWernerState(0);

fprintf("Optimizer NS2 LP");
tic
for i=1:10
p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), povms);
p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), povms);
[alpha, bellcoeffs, LPstatus, dual_alpha] = BroadcastNS2_LP(p_entangled, p_uniform);
end
toc

tic
for i=1:10
p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), povms);
p_uniform   = ProbMultidimArray(final_state(NoisyWernerState(1), channel), povms);
[alpha2, ~, ~, duals] = optimizerLP({p_entangled, p_uniform});
aux_shape = num2cell([ins, outs]);
bellcoeffs_opt = reshape(-duals{1}(1:prod([aux_shape{:}])), [aux_shape{:}]);
end
toc

fprintf("alpha=%f %f\n", alpha, alpha2);


%%

% function opt = optimizer_povm_Alice(ins, outs)
%     fprintf("Calculating YALMIP's 'optimizer' of the scenario. It can take a while.\n");
%     fprintf("Set 'OPTIMIZER_FLAG' to false to avoid 'optimizer' for debugging purposes.\n");
%     nrparties = length(ins);
% 
%     dimA = outs(1);
%     dimB1 = outs(2);
%     dimB2 = outs(3);
%     
%     final_state_dim = dimA*dimB1*dimB2;
%     final_state = sdpvar(final_state_dim,final_state_dim,'hermitian','complex');
%     ia_state = give_ia_state(final_state);
%     
%     povms_alice = cell(ins(1), outs(1));
%     partial_products = cell(ins(1), outs(1));
%     %ia_partial_products = cell(ins(1), outs(1));
%     for x=1:ins(1)
%         for a=1:outs(1)
%             povms_alice{x,a} = sdpvar(outs(1),outs(1),'hermitian','complex');
%             dimen = outs(2)*outs(3);
%             partial_products{x,a} = sdpvar(dimen,dimen,'full','complex');
%             %ia_partial_products{x,a} = partial_products{x,a}(find(triu(ones(dimen)))); 
%         end
%     end
% 
%     partyidx = 1;
%     constraints_povm = [];
%     for x=1:ins(partyidx)
%         summ = 0;
%         for a=1:outs(partyidx)
%             summ = summ + povms_alice{x,a};
%             constraints_povm = [constraints_povm,  povms_alice{x,a} >= 0];
%         end
%         constraints_povm = [constraints_povm, summ == eye(outs(partyidx))];
%     end
% 
%     
%     BellOperator = 0;
%     for x=1:ins(1)
%        for a=1:outs(1)    
%             BellOperator = BellOperator + Tensor(povms_alice{x,a}, partial_products{x,a});
%        end
%     end
%     
%    
%     objective = real(trace(final_state*BellOperator));
%     
%     opt_povm_A = optimizer(constraints_povm, -objective, ...
%                          sdpsettings('solver','mosek'), ...
%                         [{ia_state}, partial_products(:)'], ...
%                         [{objective}, povms_alice(:)']);
% 
%     opt = opt_povm_A;
% end


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