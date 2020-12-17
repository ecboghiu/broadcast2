function channels = giveLudersInstrument(dim_in, dim_out, instr_in, instr_out)

Uisometry = giveChannelRAND(dim_in, dim_out);
channels = {{{}}};
for w = 1:instr_in
    povm = RandomPOVM(dim_in, instr_out);
    kraus = cell(instr_out,1);
    for d=1:instr_out
       kraus{d} = Uisometry * sqrtm(povm{d}); 
       channels{1}{w}{d} = ChoiMatrix({kraus{d}}); 
    end
end

end

