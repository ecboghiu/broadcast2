function result = checkThatInstrChannelsAreGood(channels_w_d, instr_ins, instr_outs, dims_in, dims_out)

for instr = 1:length(instr_ins)
   for w=1:instr_ins(instr)
        channel = 0;
        for d=1:instr_outs(instr)
            channel = channel + ChoiMatrix(channels_w_d{instr}{w}{d});
        end
        dimx = dims_in(instr);
        dimy = dims_out(instr);
        tracedoutput = PartialTrace(channel, [2], [dimx, dimy]);
        diff = tracedoutput - eye(dimx);
        %assert( norm(diff) <= 1e-6, 'Choi matrix condition not satisfied.');

        choieigs = eig(channel);
        assert( all(choieigs >= -1e-6), 'Choi matrix not semi-definite positive.'); 
   end

end
result = true;
end