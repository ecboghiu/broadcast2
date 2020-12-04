function lambda_w_d = giveINSTRRandomIsometryFamily(instr_ins, instr_outs, dims_in, dims_out)
nr_instruments = length(instr_ins);
lambda_w_d = {{}};

% I generate one random isometry, and then complete up to d with
% maps which just adds a maximally mixed state for the other party
possible_placements = ["right","left"];

% TODO CHECK THAT THIS IS A VALID MAP

for instr_idx = 1:nr_instruments
    % RandomSuperoperator(DIM,TP,UN,RE,KR)
    % DIM: Either a scalar, indicating that PHI should act on DIM-by-DIM matrices, or a 1-by-2 vector, indicating that PHI should take DIM(1)-by-DIM(1) matrices as input and return DIM(2)-by-DIM(2) matrices as output.
    % TP (optional, default 1): A flag (either 1 or 0) indicating that PHI should or should not be trace-preserving.
    % UN (optional, default 0): A flag (either 1 or 0) indicating that PHI should or should not be unital (i.e., satisfy Φ(I)=I).
    % RE (optional, default 0): A flag (either 1 or 0) indicating that the Choi matrix of PHI (equivalently, its Kraus operators) should only have real entries (RE = 1) or that it is allowed to have complex entries (RE = 0).
    % KR (optional, default prod(DIM)): The maximal number of Kraus operators of the superoperator to be produced. With probability 1, it will have exactly KR Kraus operators (if KR ≤ prod(DIM)).
    
    random_weights = rand(1,instr_outs(instr_idx));
    random_weights = random_weights / sum(random_weights);
    for w=1:instr_ins(instr_idx)
        channels = {};
        channels{1} = ChoiMatrix({random_weights(1)* giveChannelRAND(dims_in(instr_idx),dims_out(instr_idx))});
        summ = channels{1};
        for d=2:instr_outs(instr_idx)
            % for everything except the first channel add the channel which
            % adds an identity to the right with some uniform weights
            % such that the channel when summed over d is unital
            placement_idx = randi([1,2]); % generate a random "right" or "left"
            channels{d} = random_weights(d)*giveChannelAddsIdentity(2,2,possible_placements(placement_idx));
            summ = summ + channels{d};
        end
        pTraceNormalization = PartialTrace(summ, 2, [dims_in(instr_idx),dims_out(instr_idx)]);
        for_assert = pTraceNormalization - diag(pTraceNormalization); 
        if norm(for_assert(:))<1e-6
           warning("Suming over outcomes doesn't give a correct Choi matrix"); 
        end
        pTraceNormalization = trace(pTraceNormalization)/dims_in(instr_idx);
        
        for d=1:instr_outs(instr_idx)
           channels{d} = channels{d}/pTraceNormalization; % now the partial trace should be normalized at the end 
           if ~IsPSD(channels{d},1e-12)
              warning("Channel not PSD");
           end
        end
        lambda_w_d{instr_idx}{w}=channels;
    end
end

end
