function out = checkPOVMsAreGood(povms, ins, outs)

for p = 1:length(ins)
   for x=1:ins(p)
       summ = 0;
      for a=1:outs(p)
          summ = summ + povms{p}{x}{a};
          assert(IsPSD(povms{p}{x}{a},1e-8),'POVM not positive!')
      end
      assert(IsPSD(summ,1e-8),"Doesn't sum to identity over outputs.");
   end
end

out = true;
end

