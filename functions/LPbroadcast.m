function funcoutput = LPbroadcast(probarray, ins, outs)  
%%%% IMPORTANT I eliminate probabilities smaller that 1e-4, THIS IS QUITE
%%%% a "big" number

party_for_det_points = 1; % this is party 'A'
    nrinputsofA  = ins(party_for_det_points);
    nroutputsofA = outs(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;
    %fprintf('nrdetpointsold = %d\n',nr_det_points);

    %alpha = sdpvar(1); % alpha will be the visibility

    %visibility_constraints = [alpha >= 0];
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
        positivityconstraints = [positivityconstraints, q(idx) >= 0];
        idx = idx + 1;
    end
    
    %fprintf('pos nr old: %d\n',length(positivityconstraints));

    % non signalling constraints
    nonsignalling_constraintsB = [];
    
    % non signaling for bob:
    for lam = 1:nr_det_points
        coordstructure = [outs(2), ins(2)];
        product = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
        for idx = 1:length(product(:,1))
            b = product(idx,1);
            y = product(idx,2);
            combinations = nchoosek(1:ins(3),2);
            for idx2 = 1:length(combinations(:,1))
                z1 = combinations(idx2,1);
                z2 = combinations(idx2,2);
                %fprintf('b y z1 z2: %d %d %d %d\n',b,y,z1,z2);
                summ1 = 0;
                for idx3 = 1:outs(3)
                    %for d = 1:outs(4)
                        summ1 = summ1 + qmatrix{lam}{y}{z1}{b}{idx3};
                    %end
                end
                
                summ2 = 0;
                for idx3 = 1:outs(3)
                    %for d = 1:outs(4)
                        summ2 = summ2 + qmatrix{lam}{y}{z2}{b}{idx3};
                    %end
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
        product = ind2subv(coordstructure, 1:prod(coordstructure,'all'));
        for idx = 1:length(product(:,1))
            c = product(idx,1);
            z = product(idx,2);
            combinations = nchoosek(1:ins(2),2);
            for idx2 = 1:length(combinations(:,1))
                y1 = combinations(idx2,1);
                y2 = combinations(idx2,2);
                summ1 = 0;
                for idx3 = 1:outs(2)
                    %for d = 1:outs(4)
                        summ1 = summ1 + qmatrix{lam}{y1}{z}{idx3}{c};
                    %end
                end
                %fprintf('c z y1 y2: %d %d %d %d\n',c,z,y1,y2);
                summ2 = 0;
                for idx3 = 1:outs(2)
                    %for d = 1:outs(4)
                        summ2 = summ2 + qmatrix{lam}{y2}{z}{idx3}{c};
                    %end
                end
                nonsignalling_constraintsC = [nonsignalling_constraintsC, summ1 == summ2];
            end
        end
    end
    %fprintf('nr NSC old: %d\n',length(nonsignalling_constraintsC));
    
    % partial normalization nonsignaling (summing over bc should not depend
    % on the inputs)
    
    nonsignalling_constraintsBC = [];
    % if I sum over all outputs this shouldn't depend on any input,
    % only on the hidden variable
    for lam = 1:nr_det_points          
        inputstructure = [ins(2) ins(3)];
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
                    %for d = 1:outs(4)
                        summ1 = summ1 + qmatrix{lam}{y1}{z1}{b}{c};
                    %end
                end
            end
            summ2 = 0;
            for b = 1:outs(2)
                for c = 1:outs(3)
                    %for d = 1:outs(4)
                        summ2 = summ2 + qmatrix{lam}{y2}{z2}{b}{c};
                    %end
                end
            end
            nonsignalling_constraintsBC = [nonsignalling_constraintsBC, summ1 == summ2];
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
          %finalstate = final_state(ini_state(alpha),channels);
          probability = probarray(x,y,z,a,b,c);%prob(finalstate,meas,[a,b,c],[x,y,z]);

          probability_constraints = [probability_constraints, summ == clean(probability,1e-6)];
          probability_constraints_inp_out = [probability_constraints_inp_out, [a;b;c;x;y;z]];
       end
    end
    
   
     
    objective = alpha;
    
    constraints = [probability_constraints:'probability', ...
                visibility_constraints:'visibility', ...
                positivityconstraints:'positivity', ...
                nonsignalling_constraintsB:'nonsignallingB',...
                nonsignalling_constraintsC:'nonsignallingC',...
                nonsignalling_constraintsBC:'nonsignallingBC'];

    optsol = optimize(constraints, -objective, ...
        sdpsettings('solver','gurobi','verbose',0,'showprogress',0, 'debug', 1, 'warning',1));

    if optsol.problem ~= 0
%         fprintf("without objective %d", optsol.problem);
%         optsol = optimize(constraints, [], ...
%         sdpsettings('solver','mosek','verbose',1,'showprogress',0, 'debug', 1, 'warning',1));

        fprintf("without objective %d", optsol.problem);
        %disp(optsol)
        break;
        %error('Check why the problem is not successfully solved.');
    else
        nrduals = length(probability_constraints);
        dualvals = zeros(nrduals,1);
        for i=1:nrduals
            dualvals(i) = value(dual(probability_constraints(i)));
        end
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

        %finalstate = final_state(ini_state(value(alpha)),channels);
        bellcoeffs(x,y,z,a,b,c) = dualvals(idx);
        %finalprob(x,y,z,a,b,c) = prob(finalstate,meas,[a,b,c],[x,y,z]);
    end
    
    funcoutput = cell(1,1);
    funcoutput{1} = bellcoeffs;
end
