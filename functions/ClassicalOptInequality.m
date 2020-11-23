function funcoutput = ClassicalOptInequality(bellcoeffs, ins, outs)
    
    party_for_det_points = 1; % this is party 'A'
    nrinputsofA  = ins(party_for_det_points);
    nroutputsofA = outs(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;

    % alpha = sdpvar(1); % alpha will be the visibility
    %finalstate = final_state(ini_state(alpha),channel);
    %visibility_constraints = [alpha >= 0, alpha <= 1];

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
    nonsignalling_constraintsB = [];
    nonsignalling_constraintsC = [];
    nonsignalling_constraintsBC = [];
    
    % non signaling for bob:
    for lam = 1:nr_det_points
        coordstructure = [outs(2) ins(2)];
        product = ind2subv(coordstructure, 1:prod(coordstructure(:)));
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
                nonsignalling_constraintsB = [nonsignalling_constraintsB,
                                            summ1 == summ2];
            end
        end
    end
    %non signaling for charlie:
    for lam = 1:nr_det_points
        coordstructure = [outs(3) ins(3)];
        product = ind2subv(coordstructure, 1:prod(coordstructure(:)));
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
                nonsignalling_constraintsC = [nonsignalling_constraintsC,
                                            summ1 == summ2];
            end
        end
    end
    
    % just in case i also add that if you summ up over b,c you get
    % something depending on lambda but not on yz
    
    %create a big vecotr with combinations of y,z
    possibleValuesForyz = zeros(ins(2)*ins(3),2);
    idx = 1;
    for y = 1:ins(2)
        for z = 1:ins(3)
            possibleValuesForyz(idx,1) = y;
            possibleValuesForyz(idx,2) = z;
            idx = idx + 1;
        end
    end
    combinations = nchoosek(1:ins(2)*ins(3),2);
    
    for lam = 1:nr_det_points
        for idx2 = 1:size(combinations,1)
            yzIDX1 = combinations(idx2,1);
            yzIDX2 = combinations(idx2,2);
            
            y1=possibleValuesForyz(yzIDX1,1);
            z1=possibleValuesForyz(yzIDX1,2);
            y2=possibleValuesForyz(yzIDX2,1);
            z2=possibleValuesForyz(yzIDX2,2);
            
            summ1 = 0;
            for b = 1:outs(2)
                for c = 1:outs(3)
                    summ1 = summ1 + qmatrix{lam}{y1}{z1}{b}{c};
                end
            end
            
            summ2 = 0;
            for b = 1:outs(2)
                for c = 1:outs(3)
                    summ1 = summ1 + qmatrix{lam}{y2}{z2}{b}{c};
                end
            end
            
            nonsignalling_constraintsBC = [nonsignalling_constraintsBC,
                                        summ1 == summ2];
        end
    end
    
    % last normalization: sum over b,c,lambda is 1 to ensure 
    % sum over abc of p(abc,xyz)=1.
    normalization_constraints = [];
    for idx = 1:ins(2)*ins(3)
        y=possibleValuesForyz(idx,1);
        z=possibleValuesForyz(idx,2);            
        summ = 0;
        for lam = 1:nr_det_points
            for b = 1:outs(2)
                for c = 1:outs(3)
                    summ = summ + qmatrix{lam}{y}{z}{b}{c};
                end
            end
        end
        normalization_constraints = [normalization_constraints, summ == 1];
    end
        
    % objective
    detA = givedetstratA(outs(1),ins(1));
    
    cartproductOUT = ind2subv(outs, 1:prod(outs(:)));
    cartproductIN  = ind2subv(ins,  1:prod(ins(:)));
    
    objective = 0;
    for i1 = 1:length(cartproductOUT)
        a = cartproductOUT(i1,1);
        b = cartproductOUT(i1,2);
        c = cartproductOUT(i1,3);
        for i2 = 1:length(cartproductIN)
            x = cartproductIN(i2,1);
            y = cartproductIN(i2,2);
            z = cartproductIN(i2,3);
            summ = 0;
            for aux = 1:nr_det_points
              summ = summ + detA(aux,x,a) * qmatrix{aux}{y}{z}{b}{c};
            end
            objective = objective + bellcoeffs(x,y,z,a,b,c)*summ;
        end
     end
  
    
    constraints = [positivityconstraints(1):'pos1', ...
                positivityconstraints(2:end):'positivity', ...
                nonsignalling_constraintsB:'nonsignallingB',...
                nonsignalling_constraintsC:'nonsignallingC',...
                nonsignalling_constraintsBC:'nonsignallingBC',...
                normalization_constraints:'normalization'];

    optsol = optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
                     'verbose',0,'dualize',0, ...
                     'showprogress',0,...
                     'savesolverinput',0,'savesolveroutput',0,'debug',0));
                
 
    funcoutput = cell(1);
    funcoutput{1} = optsol;
    funcoutput{2} = value(objective);
end

