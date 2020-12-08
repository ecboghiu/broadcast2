function probability_ndarray = ProbMultidimArrayInstrumental(inistate,povms,lambda_w_d,ins,outs)
    % povms is called as povms{party}{input}{output}
%     nrparties = max(size(povms));
%     inputs_per_party = zeros(1,nrparties);
%     outputs_per_party = zeros(1,nrparties);
%     for p=1:nrparties
%         inputs_per_party(p) = max(size(povms{p}));
%         outputs_per_party(p) = max(size(povms{p}{1}));
%     end
    inputs_per_party = ins(1:3);
    outputs_per_party = outs(1:3);
    nrparties = length(inputs_per_party);
    
    % the channel is called as lambda_w_d{instr}{w}{d}
    nrinstrs = max(size(lambda_w_d));
    instr_ins = zeros(1,nrinstrs);
    instr_outs = zeros(1,nrinstrs);
    for inst=1:nrinstrs
        instr_ins(inst) = max(size(lambda_w_d{inst}));
        instr_outs(inst) = max(size(lambda_w_d{inst}{1}));
    end
    
    
    % Define a zero multidim where to place the probability distribution
    aux = [inputs_per_party, instr_ins, outputs_per_party, instr_outs];
    dims = num2cell(aux);
    probability_ndarray = zeros(1,prod(aux));
    probability_ndarray = reshape(probability_ndarray, dims{:});
        
    
    allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
    for slice=1:size(allinputoutputcombinations,1)
        % TODO CHANGE THIS
        ins = num2cell(allinputoutputcombinations(slice,1:3)); % after nrparties we have the w
        outs = num2cell(allinputoutputcombinations(slice,5:7));

        tensor = 1;
        for p=1:nrparties
             tensor = kron(tensor,povms{p}{ins{p}}{outs{p}});
        end
        w = allinputoutputcombinations(slice,4);
        d = allinputoutputcombinations(slice,8);
        instr = 1;
        state_w_d = ApplyMapToBipartition(inistate,lambda_w_d{instr}{w}{d},"right");
        probability_ndarray(ins{:},w,outs{:},d) = real(trace(tensor*state_w_d));
    end
    %probability_ndarray = clean(probability_ndarray, 1e-6);
end
