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
state = NoisyWernerState(0.7);

optimizer_objects = optimizer_SeeSaw(dims, 2);

%%
optimizer_channel = optimizer_objects{1};
opt_a = optimizer_objects{2};
opt_b = optimizer_objects{3};
opt_c = optimizer_objects{4};

ia_povms = give_ia_povms(povms, ins, outs);
ia_state = give_ia_state(state);
ia_choi = give_ia_state(ChoiMatrix(channel));

%%
tic
for i=1:25
    [newchannel,newobjval1,~] = SeeSawOverChannel(state, bellcoeffs, povms);
end
toc
fprintf('%f\n', newobjval1);

tic
for i=1:25
    output = optimizer_channel([{ia_state}, {bellcoeffs}, ia_povms(:)']);
end
toc
fprintf('%f\n', output{2});

%%
output_state = final_state(state, channel);
partyidx = 1;
tic
for i=1:50
[newPproj,newobjval2,~] = SeeSawOverASingleParty(partyidx, output_state, bellcoeffs, povms);
end
toc
fprintf('%f\n', newobjval2);
    

output2 = opt_a([{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(2,:,:),1,[]),reshape(ia_povms(3,:,:),1,[])]);
output3 = opt_b([{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(1,:,:),1,[]),reshape(ia_povms(3,:,:),1,[])]);
output4 = opt_c([{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(1,:,:),1,[]),reshape(ia_povms(2,:,:),1,[])]);


%                         fprintf("Calculating optimizer object\n");
%     opt_choi = optimizer([constraints_choi{:}], -objective, ...
%                          sdpsettings('solver','mosek'), ...
%                         [{ia_state}, {bellcoeffs}, ia_povms(:)'], {choi, objective});

%% funcs
function opt = optimizer_SeeSaw(dims, inputdimspace)
    fprintf("Calculating YALMIP's 'optimizer' of the scenario. It can take a while.\n");
    fprintf("Set 'OPTIMIZER_FLAG' to false to avoid 'optimizer' for debugging purposes.\n");
    nrparties = length(dims)/2;
    ins = dims(1:nrparties);
    outs = dims(nrparties+1:end);

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
    
    aux_cell = num2cell(dims);
    bellcoeffs = sdpvar(aux_cell{:},'full');
    
    povms = cell(nrparties, max(ins), max(outs));
    ia_povms = cell(nrparties, max(ins), max(outs)); 
    for party=1:nrparties
        for in=1:ins(party)
            for out=1:outs(party)
                povms{party,in,out} = sdpvar(outs(party),outs(party),'hermitian','complex');
                ia_povms{party,in,out} = povms{party,in,out}(find(triu(ones(outs(party)))));
            end
        end
    end
    choi = sdpvar(choidim,choidim,'hermitian','complex');
    ia_choi = give_ia_state(choi);

    fprintf("Applying map through Choi isomorphism. \n");  
    state_small = state; %ini_state(alpha);
    state_small_reshaped = reshape(state_small, dimA,dimB,dimA,dimB);
    
   % TODO change this to a cleaner form, check that it gives the same
    state_sum = 0;
    for i=1:dimA
        for j=1:dimB
            for k=1:dimA
                for l=1:dimB
                    %scnd = PartialTrace( choi * Tensor( ketbra(j,l,dimB).', idB1B2),...
                    %                    1, [dimB,dimB1*dimB2] );
                    state_sum = state_sum + state_small_reshaped(i,j,k,l) * ...
                                    kron(ketbra(i,k,dimA), ApplyMap(ketbra(j,l,dimB), choi));
                end
            end
        end
    end
    output_state = state_sum;
  
% %     biggerstate = kron( state.', eye(dimA*dimB1*dimB2) );
% %     Phi = auxPHI(dimA);
% %     biggerchannel = kron(Phi*Phi',choi);
% %     % after previous line tensor spaces are A x A x B x B1 x B2, we need
% %     % to swap the 2nd and 3rd
% %     swapop=Tensor(eye(dimA),SwapOperator(2),eye(dimB1*dimB2)); % This ends up being a big matrix!
% %     biggerchannel = swapop * biggerchannel * swapop';
% %     output_state = PartialTrace(biggerchannel*biggerstate, [1,2], [dimA,dimB,dimA,dimB1,dimB2]);
% %     %state = AppleMap(inistate, biggerchannel);  % either this or the
% %                                                 % previous line work well
    
    fprintf("Adding constraints. \n");  
    constraints_choi = cell(2);
    constraints_choi{1} = choi >= 0;
    constraints_choi{2} = PartialTrace(choi, 2, [dimB, dimB1*dimB2]) == idB;
    
    constraints_povm = cell(nrparties, max(ins));
    for partyidx=1:nrparties
        for x=1:ins(partyidx)
            summ = 0;
            for a=1:outs(partyidx)
                summ = summ + povms{party,x,a};
                constraints_povm{partyidx,x} = [constraints_povm{partyidx,x},  povms{party,x,a} >= 0];
            end
            constraints_povm{partyidx,x} = [constraints_povm{partyidx,x}, summ == eye(outs(partyidx))];
        end
    end
    
    fprintf("Calculating Bell operator.\n");
    tic
    state_times_BellOperator = 0;
    for x=1:ins(1)
        for y=1:ins(2)
            for z=1:ins(3)
                for a=1:outs(1)
                    for b=1:outs(2)
                        for c=1:outs(3)
                            fprintf("Progress on calculating symbolically output_state*BellOperator: %d %d %d %d %d %d\n", x, y, z, a, b, c);
                            term = Tensor(povms{1,x,a}, ...
                                            povms{2,y,b},...
                                            povms{3,z,c});
                            state_times_BellOperator = state_times_BellOperator + output_state*bellcoeffs(x,y,z,a,b,c)*term;
                        end
                    end
                end
            end
        end        
    end
    toc
    
    
    
    fprintf("Calculating the objective. \n");
    objective = real(trace(state_times_BellOperator));
    

    
    fprintf("Calculating optimizer object\n");
    opt_choi = optimizer([constraints_choi{:}], -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, {bellcoeffs}, ia_povms(:)'], ...
                        {choi, objective});

    opt_povm_A = optimizer([constraints_povm{1,:,:}], -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(2,:,:),1,[]),reshape(ia_povms(3,:,:),1,[])], ...
                        [{objective}, povms(:)']);
    opt_povm_B = optimizer([constraints_povm{2,:,:}], -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(1,:,:),1,[]),reshape(ia_povms(3,:,:),1,[])], ...
                        [{objective}, povms(:)']);
    opt_povm_C = optimizer([constraints_povm{3,:,:}], -objective, ...
                         sdpsettings('solver','mosek'), ...
                        [{ia_state}, {bellcoeffs}, {ia_choi}, reshape(ia_povms(1,:,:),1,[]),reshape(ia_povms(2,:,:),1,[])], ...
                        [{objective}, povms(:)']);

    opt = {opt_choi, opt_povm_A, opt_povm_B, opt_povm_C};
                    
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