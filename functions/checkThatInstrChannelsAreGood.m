function result = checkThatInstrChannelsAreGood(channels_w_d, instr_ins, instr_outs, dims_in, dims_out)

tol = 1e-7;

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
              warning("Channels not positive! min(eig)=%g, w,d=%d,%d\n", min(eig(channels_w_d{instr}{w}{d})), w, d);
            end
        end
        tracedoutput = PartialTrace(channel, 2, [dims_in(instr), dims_out(instr)]);     
        if norm(tracedoutput-eye(dims_in(instr)),'fro')>tol
           warning("Channel not a valid Choi state! w=%d, distance to eye(d1)=%g\n", w, norm(tracedoutput-eye(dims_in(instr))));
        end
         
   end
end
result = true;
end