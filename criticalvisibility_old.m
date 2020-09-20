function funcoutput = criticalvisibility_old(meas, channels, ins, outs)  
    party_for_det_points = 1; % this is party 'A'
    nrinputsofA  = ins(party_for_det_points);
    nroutputsofA = outs(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;
    fprintf('nrdetpointsold = %d\n',nr_det_points);

    alpha = sdpvar(1); % alpha will be the visibility

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

    cartproductOUT = ind2subv(outs, 1:prod(outs,'all'));
    cartproductIN  = ind2subv(ins,  1:prod(ins,'all'));

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
    
    % nonsignaling for instrument
    % nonsignaling D not influenced by B1 and B2 (or B and C)
    nonsignalling_constraintsBC = [];
     for lam = 1:nr_det_points
        coordstructure = [outs(4) ins(4)];
        product = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
        for idx = 1:size(product,1)
            d = product(idx,1);
            w = product(idx,2);
            
            inputstructure = [ins(2), ins(3)];
            YZWinputsCartesianProduct = ind2subv(inputstructure, 1:prod(inputstructure,'all'));
            combinations = nchoosek(1:size(YZWinputsCartesianProduct,1),2);
            for index = 1:size(combinations,1)
                idx1 = combinations(index,1);
                idx2 = combinations(index,2);
                
                y1 = YZWinputsCartesianProduct(idx1,1);
                z1 = YZWinputsCartesianProduct(idx1,2);
                
                y2 = YZWinputsCartesianProduct(idx2,1);
                z2 = YZWinputsCartesianProduct(idx2,2);
                
                summ1 = 0;
                for b = 1:outs(2)
                    for c = 1:outs(3)
                        summ1 = summ1 + qmatrix{lam}{y1}{z1}{w}{b}{c}{d};
                    end
                end
                summ2 = 0;
                for b = 1:outs(2)
                    for c = 1:outs(3)
                        summ2 = summ2 + qmatrix{lam}{y2}{z2}{w}{b}{c}{d};
                    end
                end
                nonsignalling_constraintsBC = [nonsignalling_constraintsBC, ...
                                        summ1 == summ2];
            end
        end
     end

    % nonsignaling of partial normalization
    nonsignalling_constraintsBCD = [];
    % if I sum over all outputs this shouldn't depend on any input,
    % only on the hidden variable
    for lam = 1:nr_det_points          
        inputstructure = [ins(2) ins(3) ins(4)];
        YZWinputsCartesianProduct = ind2subv(inputstructure, 1:prod(inputstructure,'all'));
        combinations = nchoosek(1:size(YZWinputsCartesianProduct,1),2);
        for index = 1:size(combinations,1)
            idx1 = combinations(index,1);
            idx2 = combinations(index,2);

            y1 = YZWinputsCartesianProduct(idx1,1);
            z1 = YZWinputsCartesianProduct(idx1,2);
            w1 = YZWinputsCartesianProduct(idx1,3);

            y2 = YZWinputsCartesianProduct(idx2,1);
            z2 = YZWinputsCartesianProduct(idx2,2);
            w2 = YZWinputsCartesianProduct(idx2,3);

            summ1 = 0;
            for b = 1:outs(2)
                for c = 1:outs(3)
                    for d = 1:outs(4)
                        summ1 = summ1 + qmatrix{lam}{y1}{z1}{w1}{b}{c}{d};
                    end
                end
            end
            summ2 = 0;
            for b = 1:outs(2)
                for c = 1:outs(3)
                    for d = 1:outs(4)
                        summ2 = summ2 + qmatrix{lam}{y2}{z2}{w2}{b}{c}{d};
                    end
                end
            end
            nonsignalling_constraintsBCD = [nonsignalling_constraintsBCD, ...
                                    summ1 == summ2];
        end
    end     
     
     
     
    
     
    objective = alpha;
    
    constraints = [probability_constraints:'probability', ...
                visibility_constraints:'visibility', ...
                positivityconstraints(1):'pos1', ...
                positivityconstraints(2:end):'positivity', ...
                nonsignalling_constraintsB:'nonsignallingB',...
                nonsignalling_constraintsC:'nonsignallingC',...
                nonsignalling_constraintsBC:'nonsignallingBC',...
                nonsignalling_constraintsBCD:'nonsignallingBCD'];

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
