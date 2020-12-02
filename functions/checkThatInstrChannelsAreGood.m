function result = checkThatInstrChannelsAreGood(channels_w_d, instr_ins, instr_outs, dims_in, dims_out)

nr_instrs = size(channels_w_d,2);
nr_w = size(channels_w_d{1},2);
nr_d = size(channels_w_d{1}{1},2);
assert(nr_instrs==length(instr_ins));
assert(nr_w==instr_ins(1));
assert(nr_d==instr_outs(1));


for instr = 1:length(instr_ins)
   for w=1:instr_ins(instr)
        channel = 0;
        for d=1:instr_outs(instr)
            channel = channel + ChoiMatrix(channels_w_d{instr}{w}{d});
            assert(IsPSD(channels_w_d{instr}{w}{d},1e-8),"Not PSD");
        end
        assert(IsPSD(channel,1e-8),"Not PSD");
        tracedoutput = PartialTrace(channel, [2], [dims_in(instr), dims_out(instr)]);      
   end
end
result = true;
end