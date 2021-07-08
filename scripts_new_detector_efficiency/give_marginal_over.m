function new_prob_vec = give_marginal_over(prob_vec,sum_over,ins,outs)
% we assume the probabilities satisfy the non-signaling property
assert(max(sum_over)<=length(ins),"Shouldn't be longer than the number of parties.");

not_summed_over = setdiff(1:length(ins), sum_over);

new_prob_vec = zeros([ins(not_summed_over), outs(not_summed_over)]);

nrpartiesnotsummedover = length(not_summed_over);
inputs_per_party = ins(not_summed_over);
outputs_per_party = outs(not_summed_over);
aux = [inputs_per_party, outputs_per_party];
allinputoutputcombinations = ind2subv(aux, 1:prod(aux(:)));
for slice=1:size(allinputoutputcombinations,1)
    ins = num2cell(allinputoutputcombinations(slice,1:nrpartiesnotsummedover));
    outs = num2cell(allinputoutputcombinations(slice,nrpartiesnotsummedover+1:end));

    tensor = 1;
    for p=1:nrparties
         tensor = kron(tensor,povms{p}{ins{p}}{outs{p}});
    end

    probability_ndarray(ins{:},outs{:}) = real(trace(tensor*state));
end

end

