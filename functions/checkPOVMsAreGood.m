function out = checkPOVMsAreGood(povms, ins, outs)

tol = 1e-8;

for p = 1:length(ins)
   for x=1:ins(p)
       summ = 0;
      for a=1:outs(p)
          summ = summ + povms{p}{x}{a};
          if ~IsPSD(povms{p}{x}{a},tol)
              warning("POVM not positive! min eig %g", min(eig(povms{p}{x}{a})));
          end
      end
      if ~IsPSD(summ,tol)
          warning("POVM not positive! min eig %g", min(eig(povms{p}{x}{a})));
      end
   end
end

out = true;
end

