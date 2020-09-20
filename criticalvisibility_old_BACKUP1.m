function funcoutput = criticalvisibility_old(meas, channels, ins, outs) 
    % ins is something like {[1,2,3],[1,2],[1,2]} giving the 
    % inputs of each party. cell 1 is for party A, cell 2 is for party B
    % and cell 3 is for party C. Similarly outs is something of the
    % form {[1,2],[1,2],[1,2,3,4]} each cell element giving the different
    % outputs of each party.
    % in practice we only care about the NUMBER of inputs or outputs
    % and never about what they are exactly
    
    
    party_for_det_points = 1; % this is party 'A'
    nrinputsofA  = ins(party_for_det_points);
    nroutputsofA = outs(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;
    fprintf('nrdetpointsold = %d\n',nr_det_points);
 

    alpha = sdpvar(1); % alpha will be the visibility

    visibility_constraints = [alpha >= 0, alpha <= 1];


    q = sdpvar(nr_det_points * prod([ins,outs])/nrinputsofA/nroutputsofA, 1);    

    
    % alternative to what follows: define an index function which to 
    % every element i of array q assigns the element (lam,y,z,b,c) where
    % we interpret (lam,y,z,b,c) is a decomposition of i according to
    % i = c+b*2+z*2*2+y*2*2*2+lam*2*2*2*2
    % this can be easily achieved with something like reshape but I chose
    % cells for this
    positivityconstraints = [];
%     qmatrix = {{{{{{{{0}}}}}}}};
%     idx = 1;
%     for lam = 1:nr_det_points
%        for y = ins{2}
%            for z = ins{3}
%                for w = ins{4}
%                    for b = outs{2}
%                    for c = outs{3}
%                        for d = outs{4}
%                        qmatrix{lam}{y}{z}{w}{b}{c}{d} = q(idx);
%                        positivityconstraints = [positivityconstraints, ...
%                                                 q(idx) >= 0, q(idx) <= 1];
%                        idx = idx + 1;
%                        end
%                    end
%                    end
%                end
%            end
%        end
%     end
varnumbers = [nr_det_points, ins(2:end), outs(2:end)];
    loopvars = ind2subv(varnumbers, 1:prod(varnumbers,'all'));
    
    idx = 1;
    for i=1:size(loopvars,1)
        lam = loopvars(i,1);
        y   = loopvars(i,2);
        z   = loopvars(i,3);
        w   = loopvars(i,4);
        b   = loopvars(i,5);
        c   = loopvars(i,6);
        d   = loopvars(i,7);

        qmatrix{lam}{y}{z}{w}{b}{c}{d} = q(idx);
        positivityconstraints = [positivityconstraints, ...
                                q(idx) >= 0, q(idx) <= 1];
        
        idx = idx + 1;
    end
    


    fprintf('pos nr old: %d\n',length(positivityconstraints));
    
    %normalizationconstraint = [q(1)>=0];
%     summ = 0;
%     for lam = 1:nr_det_points
%        for y = ins{2}
%            for z = ins{3}
%                for b = outs{2}
%                    for c = outs{3}
%                        summ = summ + qmatrix{lam}{y}{z}{b}{c};
%                    end
%                end
%            end
%        end
%     end
%     normalizationconstraint = [summ == 1];


    % non signalling constraints
    nonsignalling_constraintsB = [];
    
    % non signaling for bob:
    for lam = 1:nr_det_points
        coordstructure = [outs(2), ins(2), ins(4)];
        product = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
        for idx = 1:length(product(:,1))
            b = product(idx,1);
            y = product(idx,2);
            w = product(idx,3);
            combinations = nchoosek(1:ins(3),2);
            for idx2 = 1:length(combinations(:,1))
                z1 = combinations(idx2,1);
                z2 = combinations(idx2,2);
                %fprintf('b y z1 z2: %d %d %d %d\n',b,y,z1,z2);
                summ1 = 0;
                for idx3 = 1:outs(3)
                    for d = 1:outs(4)
                    summ1 = summ1 + qmatrix{lam}{y}{z1}{w}{b}{idx3}{d};
                    end
                end
                
                summ2 = 0;
                for idx3 = 1:outs(3)
                    for d = 1:outs(4)
                    summ2 = summ2 + qmatrix{lam}{y}{z2}{w}{b}{idx3}{d};
                    end
                end
                nonsignalling_constraintsB = [nonsignalling_constraintsB,
                                            summ1 == summ2];
            end
        end
    end
    fprintf('nr NSB old: %d\n',length(nonsignalling_constraintsB));

    %non signaling for charlie:
    nonsignalling_constraintsC = [];
    for lam = 1:nr_det_points
        coordstructure = [outs(3), ins(3), ins(4)];
        product = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
        for idx = 1:length(product(:,1))
            c = product(idx,1);
            z = product(idx,2);
            w = product(idx,3);
            combinations = nchoosek(1:ins(2),2);
            for idx2 = 1:length(combinations(:,1))
                y1 = combinations(idx2,1);
                y2 = combinations(idx2,2);
                summ1 = 0;
                for idx3 = 1:outs(2)
                    for d = 1:outs(4)
                    summ1 = summ1 + qmatrix{lam}{y1}{z}{w}{idx3}{c}{d};
                    end
                end
                %fprintf('c z y1 y2: %d %d %d %d\n',c,z,y1,y2);
                summ2 = 0;
                for idx3 = 1:outs(2)
                    for d = 1:outs(4)
                    summ2 = summ2 + qmatrix{lam}{y2}{z}{w}{idx3}{c}{d};
                    end
                end
                nonsignalling_constraintsC = [nonsignalling_constraintsC,
                                        summ1 == summ2];
            end
        end
    end
    fprintf('nr NSC old: %d\n',length(nonsignalling_constraintsC));

    
    
    
    % probability constraints
    det = givedetstratA(outs(1),ins(1));

    probability_constraints = [];
    probability_constraints_inp_out = [];
%                  ins = {[1,2,3],[1,2],[1,2],[1]};
% outs = {[1,2],[1,2],[1,2],[1]};   
    outs2=[2,2,2,1];
    ins2 = [3,2,2,1];
    coordstructure = [outs2(1), outs2(2), outs2(3), outs2(4)];
    cartproductOUT = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
coordstructure = [ins2(1), ins2(2), ins2(3), ins2(4)];
    cartproductIN = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
    for i1 = 1:length(cartproductOUT)
       for i2 = 1:length(cartproductIN)
          a = cartproductOUT(i1,1);
          b = cartproductOUT(i1,2);
          c = cartproductOUT(i1,3);
          d = cartproductOUT(i1,4);
          x = cartproductIN(i2,1);
          y = cartproductIN(i2,2);
          z = cartproductIN(i2,3);
          w = cartproductIN(i2,4);
          %fprintf('a b c x y z: %d %d %d %d %d %d\n',a,b,c,x,y,z);

          summ = 0;
          for aux = 1:nr_det_points
              summ = summ + det(aux,x,a) * qmatrix{aux}{y}{z}{w}{b}{c}{d};
          end
          finalstate = final_state(ini_state(alpha),channels{w}{d});
          probability = prob(finalstate,meas,[a,b,c],[x,y,z]);
          probability_constraints = [probability_constraints, ...
                                        summ == probability];
          probability_constraints_inp_out = ...
              [probability_constraints_inp_out, [a;b;c;d;x;y;z;w]];
       end
    end

    objective = alpha;
    
    constraints = [probability_constraints:'probability', ...
                visibility_constraints:'visibility', ...
                positivityconstraints(1):'pos1', ...
                positivityconstraints(2:end):'positivity', ...
                nonsignalling_constraintsB:'nonsignallingB',...
                nonsignalling_constraintsC:'nonsignallingC'];

    optsol = optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
    'verbose',0,'dualize',0, ...
    'showprogress',0,...
    'savesolverinput',0,'savesolveroutput',0,'debug',0));


    nrduals = length(probability_constraints);
    dualvals = zeros(nrduals,1);
    for i=1:nrduals
        dualvals(i) = value(dual(probability_constraints(i)));
    end
    
%     bellcoeffs = zeros(max(ins{1}),max(ins{2}),max(ins{3}),max(ins{4}), ...
%                     max(outs{1}),max(outs{2}),max(outs{3}),max(outs{4}));
%                 
%     finalprob = zeros(max(ins{1}),max(ins{2}),max(ins{3}),max(ins{4}), ...
%                     max(outs{1}),max(outs{2}),max(outs{3}),max(outs{4}));
    dimcell = num2cell([ins outs]);
    bellcoeffs = zeros(dimcell{:});           
    finalprob  = zeros(dimcell{:});
    finalAlpha = value(alpha);
    for idx=1:length(probability_constraints_inp_out(1,:))
        vec = probability_constraints_inp_out(:,idx);
        a = vec(1);
        b = vec(2);
        c = vec(3);
        d = vec(4);
        x = vec(5);
        y = vec(6);
        z = vec(7);
        w = vec(8);

        finalstate = final_state(ini_state(value(alpha)),channels{w}{d});
        bellcoeffs(x,y,z,w,a,b,c,d) = dualvals(idx);
        finalprob(x,y,z,w,a,b,c,d) = prob(finalstate,meas,[a,b,c],[x,y,z]);
    end
    
    funcoutput = cell(1);
    funcoutput{1} = finalAlpha;
    funcoutput{2} = bellcoeffs;
    funcoutput{3} = finalprob;
    funcoutput{4} = constraints;
end
