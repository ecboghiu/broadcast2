ins=[3,2,2];
outs=[2,2,2];

load('InfeasibleLP_channels.mat');
load('InfeasibleLP_meas.mat');



party_for_det_points = 1; % this is party 'A'
nrinputsofA  = ins(party_for_det_points);
nroutputsofA = outs(party_for_det_points);

nr_det_points = nroutputsofA^nrinputsofA;

alpha = sdpvar(1); % alpha will be the visibility
alpha=1;

probarray = clean(giveProbNDArray(final_state(ini_state(alpha),channels),meas,ins,outs),1e-6);

visibility_constraints = [alpha >= 0, alpha <= 1];
q = sdpvar(nr_det_points * prod([ins,outs])/nrinputsofA/nroutputsofA, 1);    

positivityconstraints = [];
varnumbers = [nr_det_points, ins(2:end), outs(2:end)];
loopvars = ind2subv(varnumbers, 1:prod(varnumbers,'all'));
idx = 1;

for i=1:size(loopvars,1)
    lam = loopvars(i,1);
    y   = loopvars(i,2);
    z   = loopvars(i,3);
    b   = loopvars(i,4);
    c   = loopvars(i,5);

    qmatrix{lam}{y}{z}{b}{c} = q(idx);
    positivityconstraints = [positivityconstraints, ...
                            q(idx) >= 0];
    idx = idx + 1;
end

% non signalling constraints
nonsignalling_constraintsB = [];

% non signaling for bob:
for lam = 1:nr_det_points
    coordstructure = [outs(2), ins(2)];
    producto = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
    for idx = 1:length(producto(:,1))
        b = producto(idx,1);
        y = producto(idx,2);
        combinaciones = nchoosek(1:ins(3),2);
        for idx2 = 1:size(combinaciones,1)
            z1 = combinaciones(idx2,1);
            z2 = combinaciones(idx2,2);
            summ1 = 0;
            for c = 1:outs(3)
                summ1 = summ1 + qmatrix{lam}{y}{z1}{b}{c};
            end
            summ2 = 0;
            for c = 1:outs(3)
                summ2 = summ2 + qmatrix{lam}{y}{z2}{b}{c};
            end
            nonsignalling_constraintsB = [nonsignalling_constraintsB, summ1 == summ2];
        end
    end
end
%fprintf('nr NSB old: %d\n',length(nonsignalling_constraintsB));

%non signaling for charlie:
nonsignalling_constraintsC = [];
for lam = 1:nr_det_points
    coordstructure = [outs(3), ins(3)];
    producto = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
    for idx = 1:length(producto(:,1))
        c = producto(idx,1);
        z = producto(idx,2);
        combinaciones = nchoosek(1:ins(2),2);
        for idx2 = 1:size(combinaciones,1)
            y1 = combinaciones(idx2,1);
            y2 = combinaciones(idx2,2);
            summ1 = 0;
            for b = 1:outs(2)
                summ1 = summ1 + qmatrix{lam}{y1}{z}{b}{c};
            end
            summ2 = 0;
            for b = 1:outs(2)
                summ2 = summ2 + qmatrix{lam}{y2}{z}{b}{c};
            end
            nonsignalling_constraintsC = [nonsignalling_constraintsC, summ1 == summ2];
        end
    end
end

% probability constraints
det = givedetstratA(outs(1),ins(1));

probability_constraints = [];
probability_constraints_inp_out = [];

cartproductOUT = ind2subv(outs, 1:prod(outs,'all'));
cartproductIN  = ind2subv(ins,  1:prod(ins,'all'));

for i1 = 1:length(cartproductOUT)
   for i2 = 1:length(cartproductIN)
      a = cartproductOUT(i1,1);
      b = cartproductOUT(i1,2);
      c = cartproductOUT(i1,3);
      x = cartproductIN(i2,1);
      y = cartproductIN(i2,2);
      z = cartproductIN(i2,3);
      %fprintf('a b c x y z: %d %d %d %d %d %d\n',a,b,c,x,y,z);

      summ = 0;
      for aux = 1:nr_det_points
          summ = summ + det(aux,x,a) * qmatrix{aux}{y}{z}{b}{c};
      end
      finalstate = final_state(ini_state(alpha),channels);
      probability = prob(finalstate,meas,[a,b,c],[x,y,z]);
      probability = probarray(x,y,z,a,b,c);
      probability_constraints = [probability_constraints, ...
                                    summ == clean(probability,1e-8)];
      probability_constraints_inp_out = ...
          [probability_constraints_inp_out, [a;b;c;x;y;z]];
   end
end



objective = alpha;

constraints = [probability_constraints:'probability', ...
            positivityconstraints(1):'positivity', ...
            nonsignalling_constraintsB:'nonsignallingB',...
            nonsignalling_constraintsC:'nonsignallingC'];

optsol = optimize(constraints,-objective, ...
    sdpsettings('solver','gurobi'));%,'verbose',1,'showprogress',1, 'debug', 1, 'warning',1));

% if optsol.problem ~= 0
%     fprintf("with objective %d", optsol.problem);
%     optsol = optimize(constraints, [], ...
%     sdpsettings('solver','gurobi','verbose',1,'showprogress',0, 'debug', 1, 'warning',1));
% 
%     fprintf("without objective %d", optsol.problem);
%     disp(optsol)
%     error('Check why the problem is not successfully solved.');
% end

nrduals = length(probability_constraints);
dualvals = zeros(nrduals,1);
for i=1:nrduals
    dualvals(i) = value(dual(probability_constraints(i)));
end

dimcell = num2cell([ins outs]);
bellcoeffs = zeros(dimcell{:});           
finalprob  = zeros(dimcell{:});
finalAlpha = value(alpha);
for idx=1:length(probability_constraints_inp_out(1,:))
    vec = probability_constraints_inp_out(:,idx);
    a = vec(1);
    b = vec(2);
    c = vec(3);
    x = vec(4);
    y = vec(5);
    z = vec(6);

    finalstate = final_state(ini_state(value(alpha)),channels);
    bellcoeffs(x,y,z,a,b,c) = dualvals(idx);
    finalprob(x,y,z,a,b,c) = prob(finalstate,meas,[a,b,c],[x,y,z]);
end

funcoutput = cell(1,4);
funcoutput{1} = finalAlpha;
funcoutput{2} = bellcoeffs;
funcoutput{3} = finalprob;
funcoutput{4} = constraints;