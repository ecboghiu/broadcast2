function lambda_w_d = giveINSTRRandomIsometryFamily(instr_ins, instr_outs, dims_in, dims_out)
nr_instruments = length(instr_ins);
lambda_w_d = {{}};
possible_placements = ["right","left"];

% TODO CHECK THAT THIS IS A VALID MAP

for instr_idx = 1:nr_instruments
    % RandomSuperoperator(DIM,TP,UN,RE,KR)
    % DIM: Either a scalar, indicating that PHI should act on DIM-by-DIM matrices, or a 1-by-2 vector, indicating that PHI should take DIM(1)-by-DIM(1) matrices as input and return DIM(2)-by-DIM(2) matrices as output.
    % TP (optional, default 1): A flag (either 1 or 0) indicating that PHI should or should not be trace-preserving.
    % UN (optional, default 0): A flag (either 1 or 0) indicating that PHI should or should not be unital (i.e., satisfy Φ(I)=I).
    % RE (optional, default 0): A flag (either 1 or 0) indicating that the Choi matrix of PHI (equivalently, its Kraus operators) should only have real entries (RE = 1) or that it is allowed to have complex entries (RE = 0).
    % KR (optional, default prod(DIM)): The maximal number of Kraus operators of the superoperator to be produced. With probability 1, it will have exactly KR Kraus operators (if KR ≤ prod(DIM)).
    channels = {};
    random_weights = rand(1,instr_outs(instr_idx));
    random_weights = random_weights / sum(random_weights);
    for w=1:instr_ins(instr_idx)
        channels{1} = ChoiMatrix({random_weights(1) * giveChannelRAND(dims_in(instr_idx),dims_out(instr_idx))});
        for d=2:instr_outs(instr_idx)
            % for everything except the first channel add the channel which
            % adds an identity to the right with some uniform weights
            % such that the channel when summed over d is unital
            placement_idx = randi([1,2]); % generate a random "right" or "left"
            channels{d} = giveChannelAddsIdentity(2,2,possible_placements(placement_idx),random_weights(d));
        end
        lambda_w_d{instr_idx}{w}=channels;
    end
end

end
