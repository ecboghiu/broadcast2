mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
newdir3 = strcat(newdir,filesep,'scripts',filesep);
addpath(newdir);
addpath(newdir2);
addpath(newdir3);

load('best_for_sqrt3.mat');
%[newalpha, newbellcoeffs, LPstatus, dual_alpha] = BroadcastLP(p_entangled, p_uniform);

channel = best_channels{1};
POVMs = best_povm{1};

ins = [3,2,2];
outs = [2,2,2];
nr_parties = length(ins);
efficiencies = [1,1,1];

p_entangled = ProbMultidimArray(final_state(NoisyWernerState(0), channel), POVMs);

%% Check normalization
efficiencies = rand(1,3);
newp = giveDetectorEfficiencydistrib(p_entangled, efficiencies, ins, outs);
outs_plus_one = outs + 1;
inputs_combs = ind2subv(ins, 1:prod(ins(:)));
outputs_combs = ind2subv(outs_plus_one, 1:prod(outs_plus_one(:)));
for input_slice=1:size(inputs_combs,1)
    inputs = num2cell(inputs_combs(input_slice,:));
    summ = 0;
    for output_slice=1:size(outputs_combs,1)
        outputs = num2cell(outputs_combs(output_slice,:));
        summ = summ + newp(inputs{:}, outputs{:});
    end
    assert(abs(summ-1)<1e-6, "Not normalized");
end

%% For efficiency 1
efficiencies = [1,1,1];
newp = giveDetectorEfficiencydistrib(p_entangled, efficiencies, ins, outs);

aux = [ins, outs+1];
allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
for slice=1:size(allinputoutputcombinations,1)
    ins_slice = num2cell(allinputoutputcombinations(slice,1:nr_parties));
    outs_slice = num2cell(allinputoutputcombinations(slice,nr_parties+1:end));

    outs_num = [outs_slice{:}];
    if any(outs_num == outs+1)
        assert(newp(ins_slice{:},outs_slice{:}) == 0, "Should be zero for perfect efficiency");
    end
end

%% For other efficiencies
% We know that it should be 
% p(abc|xyz) = eta_A p(abc|xyz) + (1-eta_A) p(bc|yz)

efficiencies = [0.9,0.8,0.7];
newp = giveDetectorEfficiencydistrib(p_entangled, efficiencies, ins, outs);

aux = [ins, outs+1];
allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
for slice=1:size(allinputoutputcombinations,1)
    ins_slice = num2cell(allinputoutputcombinations(slice,1:nr_parties));
    outs_slice = num2cell(allinputoutputcombinations(slice,nr_parties+1:end));

    outs_num = [outs_slice{:}];
    if outs_num(1) == outs(1)+1 && outs_num(2) ~= outs(2)+1 && outs_num(3) ~= outs(3)+1
        marginal = sum(p_entangled, 3 + 1);
        new_outs_slice = outs_slice;
        new_outs_slice{1} = 1;
        assert(abs(    newp(ins_slice{:},outs_slice{:}) - ...
                    (1-efficiencies(1))*efficiencies(2)*efficiencies(3)*marginal(ins_slice{:}, new_outs_slice{:}) ...
                    ) < 1e-6);
    end
    
    if outs_num(1) ~= outs(1)+1 && outs_num(2) == outs(2)+1 && outs_num(3) ~= outs(3)+1
        marginal = sum(p_entangled, 3 + 2);
        new_outs_slice = outs_slice;
        new_outs_slice{2} = 1;
        assert(abs(newp(ins_slice{:},outs_slice{:}) - ...
                    (1-efficiencies(2))*efficiencies(1)*efficiencies(3)*marginal(ins_slice{:}, new_outs_slice{:})) < 1e-6, "Should be zero for perfect efficiency");
    end
    
    if outs_num(1) ~= outs(1)+1 && outs_num(2) ~= outs(2)+1 && outs_num(3) == outs(3)+1
        marginal = sum(p_entangled, 3 + 3);
        new_outs_slice = outs_slice;
        new_outs_slice{3} = 1;
        assert(abs(newp(ins_slice{:},outs_slice{:}) - ...
                    (1-efficiencies(3))*efficiencies(1)*efficiencies(2)*marginal(ins_slice{:}, new_outs_slice{:})) < 1e-6, "Should be zero for perfect efficiency");
    end
    
    if outs_num(1) ~= outs(1)+1 && outs_num(2) == outs(2)+1 && outs_num(3) == outs(3)+1
        marginal = sum(p_entangled, [3 + 2, 3 + 3]);
        new_outs_slice = outs_slice;
        new_outs_slice{2} = 1;
        new_outs_slice{3} = 1;
        assert(abs(newp(ins_slice{:},outs_slice{:}) - ...
                    (1-efficiencies(2))*(1-efficiencies(3))*efficiencies(1)*marginal(ins_slice{:}, new_outs_slice{:})) < 1e-6, "Should be zero for perfect efficiency");
    end
end