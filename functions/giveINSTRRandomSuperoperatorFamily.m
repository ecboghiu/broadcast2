function lambda_w_d = giveINSTRRandomSuperoperatorFamily(instr_ins, instr_outs, dims_in, dims_out)
nr_instruments = length(instr_ins);
lambda_w_d = {{}};
for instr_idx = 1:nr_instruments
    % RandomSuperoperator(DIM,TP,UN,RE,KR)
    % DIM: Either a scalar, indicating that PHI should act on DIM-by-DIM matrices, or a 1-by-2 vector, indicating that PHI should take DIM(1)-by-DIM(1) matrices as input and return DIM(2)-by-DIM(2) matrices as output.
    % TP (optional, default 1): A flag (either 1 or 0) indicating that PHI should or should not be trace-preserving.
    % UN (optional, default 0): A flag (either 1 or 0) indicating that PHI should or should not be unital (i.e., satisfy Φ(I)=I).
    % RE (optional, default 0): A flag (either 1 or 0) indicating that the Choi matrix of PHI (equivalently, its Kraus operators) should only have real entries (RE = 1) or that it is allowed to have complex entries (RE = 0).
    % KR (optional, default prod(DIM)): The maximal number of Kraus operators of the superoperator to be produced. With probability 1, it will have exactly KR Kraus operators (if KR ≤ prod(DIM)).
    dims = [dims_in(instr_idx),dims_out(instr_idx)];
    for w=1:instr_ins(instr_idx)
        d = instr_outs(instr_idx);
        kraus = KrausOperators(RandomSuperoperator(dims,1,0,0,d),dims);
        ToSmallChoi = {};
        for kr=1:length(kraus)
            ToSmallChoi{kr} = ChoiMatrix({kraus{kr}});
        end
        lambda_w_d{instr_idx}{w} = ToSmallChoi;
    end
end

end