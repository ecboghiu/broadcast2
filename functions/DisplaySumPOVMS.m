function [outputArg1,outputArg2] = DisplaySumPOVMS(Pproj, ins, outs)
    for partyidx = [1,2,3]
        for x=ins{partyidx}
            fprintf("result for p=%d x=%d\n", partyidx, x);
            disp([newPproj{partyidx}{x}{1},newPproj{partyidx}{x}{2},newPproj{partyidx}{x}{1}+newPproj{partyidx}{x}{2}]); 
        end
    end
end

