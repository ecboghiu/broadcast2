function out = LocalBoundOfIneq(bellcoeffs,ins,outs)

% just for compatibility with older code with different inputs
if isa(ins,'cell') && isa(outs,'cell')
   ins = [max(ins{1}),max(ins{2}),max(ins{3})];
   outs = [max(outs{1}),max(outs{2}),max(outs{3})];
end

% This is in the broadcasting scenario of 
% oseph Bowles, Flavien Hirsch, and Daniel Cavalcanti.
% Single-copy activation of bellnonlocality via broadcasting of quantum states.arXivpreprintarXiv:2007.16034, 2020
% Therefore we have 3 parties A and B1 and B2. The local model is
% p(ab1b2|xy1y2) = \sum_l D(a|x lam) q_{b1 b2 y1 y2 lam}
% with more constraints.

party_for_det_points = 1; % this is party 'A'

nrinputsofA  = ins(party_for_det_points);
nroutputsofA = outs(party_for_det_points);

nr_det_points = nroutputsofA^nrinputsofA;

q = sdpvar(nr_det_points*ins(2)*ins(3)*outs(2)*outs(3),1);


% alternative to what follows: define an index function which to 
% every element i of array q assigns the element (lam,y,z,b,c) where
% we interpret (lam,y,z,b,c) is a decomposition of i according to
% i = c+b*2+z*2*2+y*2*2*2+lam*2*2*2*2
% this can be easily achieved with something like reshape but I chose
% cells for this
positivityconstraints = [];
qmatrix = {{{{{{0}}}}}};
idx = 1;
for lam = 1:nr_det_points
   for y = 1:ins(2)
       for z = 1:ins(3)
           for b = 1:outs(2)
               for c = 1:outs(3)
                   qmatrix{lam}{y}{z}{b}{c} = q(idx);
                   positivityconstraints = [positivityconstraints, ...
                                            q(idx) >= 0];
                   idx = idx + 1;
               end
           end
       end
   end
end

% non signalling constraints
nonsignalling_constraints = [];

% non signaling for bob:
for lam = 1:nr_det_points
    coordstructure = [outs(2) ins(2)];
    product = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
    for idx = 1:size(product,1)
        b = product(idx,1);
        y = product(idx,2);
        combinations = nchoosek(1:ins(3),2);
        for idx2 = 1:size(combinations,1)
            z1 = combinations(idx2,1);
            z2 = combinations(idx2,2);
            summ1 = 0;
            for idx3 = 1:outs(3)
                summ1 = summ1 + qmatrix{lam}{y}{z1}{b}{idx3};
            end

            summ2 = 0;
            for idx3 = 1:outs(3)
                summ2 = summ2 + qmatrix{lam}{y}{z2}{b}{idx3};
            end
            nonsignalling_constraints = [nonsignalling_constraints,
                                        summ1 == summ2];
        end
    end
end
%non signaling for charlie:
for lam = 1:nr_det_points
    coordstructure = [outs(3) ins(3)];
    product = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
    for idx = 1:size(product,1)
        c = product(idx,1);
        z = product(idx,2);
        combinations = nchoosek(1:ins(2),2);
        for idx2 = 1:size(combinations,1)
            y1 = combinations(idx2,1);
            y2 = combinations(idx2,2);
            summ1 = 0;
            for idx3 = 1:outs(2)
                summ1 = summ1 + qmatrix{lam}{y1}{z}{idx3}{c};
            end

            summ2 = 0;
            for idx3 = 1:outs(2)
                summ2 = summ2 + qmatrix{lam}{y2}{z}{idx3}{c};
            end
            nonsignalling_constraints = [nonsignalling_constraints,
                                        summ1 == summ2];
        end
    end
end

% normalization constraints for b1 b2
normb1b2 = [];
elslist = {};
iter = 1;
for lam = 1:nr_det_points
    for y = 1:ins(2)
        for z = y:ins(3)
            summ = 0;
            for b = 1:outs(2)
            for c = 1:outs(3)
                summ = summ + qmatrix{lam}{y}{z}{b}{c};
            end
            end
            
            if iter == 1
                elslist{1} = summ;
            else
                normb1b2 = [normb1b2, summ == elslist{iter-1}];
            end
            
        end
    end
end


% probability constraints
det = givedetstratA(1:outs(1),1:ins(1),nr_det_points);

dims = num2cell([ins,outs]);
prob = zeros(dims{:});

cartproductOUT = ind2subv(outs, 1:prod(outs,'all'));
cartproductIN  = ind2subv(ins,  1:prod(ins,'all'));

objective = 0;
for i1 = 1:size(cartproductOUT,1)
   a = cartproductOUT(i1,1);
   b = cartproductOUT(i1,2);
   c = cartproductOUT(i1,3);
   for i2 = 1:size(cartproductIN,1)
      x = cartproductIN(i2,1);
      y = cartproductIN(i2,2);
      z = cartproductIN(i2,3);
      summ = 0;
      for aux = 1:nr_det_points
          summ = summ + det(aux,x,a) * qmatrix{aux}{y}{z}{b}{c};
      end
      prob(x,y,z,a,b,c) = summ;
      objective = objective + prob(x,y,z,a,b,c) * bellcoeffs(x,y,z,a,b,c);
   end
end

% just in case even more normalization constraints
for i1 = 1:size(cartproductOUT,1)
   a = cartproductOUT(i1,1);
   b = cartproductOUT(i1,2);
   c = cartproductOUT(i1,3);
   for i2 = 1:size(cartproductIN,1)
      x = cartproductIN(i2,1);
      y = cartproductIN(i2,2);
      z = cartproductIN(i2,3);
      summ = 0;
      for aux = 1:nr_det_points
          summ = summ + prob(x,y,z,a,b,c);
      end
   end
end
finalnormconstraints = [summ == 1];



constraints = [positivityconstraints:'positivity',e...
                nonsignalling_constraints:'nonsignalling',...
                finalnormconstraints:'finalnorm',...
                normb1b2:'iterb1b2'];

optsol = optimize(constraints, -objective, ...
    sdpsettings('solver','mosek', ...
'verbose',0,'dualize',0, ...
'showprogress',0,...
'savesolverinput',0,'savesolveroutput',0,'debug',0));

out = cell(1);
out{1} = value(objective);
out{2} = optsol;
end

