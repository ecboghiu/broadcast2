function channels = giveLudersInstrument(dim_in, dim_out, instr_in, instr_out)

Uisometry = giveChannelRAND(dim_in, dim_out);
channels = {{{}}};
for w = 1:instr_in
    povm = RandomPOVM(dim_in, instr_out);
    for d=1:instr_out
       channels{1}{w}{d} = ChoiMatrix({Uisometry * sqrt(povm{d})}); 
    end
    %channel = ChoiMatrix(kraus);
    %channels{1}{w} = channel;
end


end

