function out = checkPOVMsAreGood(povms, ins, outs)

for p = 1:length(ins)
   for x=1:ins(p)
       summ = 0;
      for a=1:outs(p)
          summ = summ + povms{p}{x}{a};
          eigvals=eig(povms{p}{x}{a});
          assert(all(eigvals >= -1e-7),'POVM not positive!')
      end
      diff = summ-eye(outs(p));
      assert(norm(diff)<1e-7,"Doesn't sum to identity over outputs.");
      
   end
end

out = true;
end

