function result = checkThatInstrChannelsAreGood(channels_w_d, instr_ins, instr_outs, dims_in, dims_out)

tol = 1e-8;

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
            if ~IsPSD(channels_w_d{instr}{w}{d},tol)
              warning("POVM not positive! min eig %g", min(eig(channels_w_d{instr}{w}{d})));
            end
        end
        if ~IsPSD(channel,tol)
           warning("POVM not positive! min eig %g", min(eig(channel)));
        end
        tracedoutput = PartialTrace(channel, [2], [dims_in(instr), dims_out(instr)]);      
   end
end
result = true;
end