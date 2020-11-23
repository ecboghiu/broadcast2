function marginal = GiveMarginals(prob, ins, outs, idx_over_which_to_marginalize)

sizes=size(prob);
nrparties=sizes/2;

marginal = {{{{{}}}}};


for i1 = 1:length(cartproductOUT)
   for i2 = 1:length(cartproductIN)
      a = cartproductOUT(i1,1);
      b = cartproductOUT(i1,2);
      c = cartproductOUT(i1,3);
      x = cartproductIN(i2,1);
      y = cartproductIN(i2,2);
      z = cartproductIN(i2,3);

      summ = 0;
      for aux = 1:nr_det_points
          summ = summ + det(aux,x,a) * qmatrix{aux}{y}{z}{b}{c};
      end

      probability = prob(finalstate,meas,[a,b,c],[x,y,z]);
      probability_constraints = [probability_constraints, ...
                                    summ == probability];
      probability_constraints_inp_out = ...
          [probability_constraints_inp_out, [a;b;c;x;y;z]];
   end
end

end

