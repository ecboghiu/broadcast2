function [final_alpha,bellcoeffs,LPstatus,dual_alpha] = BroadcastInstrumentLP(p1,p2,nr_inputs_per_party,nr_outputs_per_party)
    % IMPORTANT alpha is defined as (1-alpha)*p1 + (alpha)*p2
    % usually p1 is the entangled, and p2 the uniform, so this
    % program minimizes alpha using this convention
    % if p1 and p2 are both local then all alpha are good
    % if p1 and p2 are both nonlocal the program is infeasible
    % ideally p1 should be nonlocal

    % IMPORTANT I use the convention where the probabilities are called as
    % p(x,y,z,a,b,c) instead of p(x,y,z,w,a,b,c,d)
    dims_p = size(p1);
    
%    assert(mod(length(dims_p),2) == 0, "The probability array should have equal number of inputs and outputs.");
    assert(all(size(p1) == size(p2)), "The two probability arrays should have equal dimenisons.");
    if norm( p1(:) - p2(:) ) < 1e-6
        disp("The probabilities should not be almost equal.");
        LPstatus = 1;
        aux_cell = num2cell(dims_p);
        bellcoeffs = zeros(aux_cell{:});
        dual_alpha=0;
        final_alpha = 0;
        return;
    end
    
    nrparties = length(nr_inputs_per_party);
    
    %nr_inputs_per_party = dims_p(1:nrparties);
    %nr_outputs_per_party = dims_p(nrparties+1:end);

    party_for_det_points = 1; % this is party 'A' % TODO Make this a function input
    nrinputsofA  = nr_inputs_per_party(party_for_det_points);
    nroutputsofA = nr_outputs_per_party(party_for_det_points);
    
    nr_det_points = nroutputsofA^nrinputsofA;
  
    alpha = sdpvar(1); % alpha will be the visibility

    tempdims = [nr_det_points, nr_inputs_per_party(2:end), nr_outputs_per_party(2:end)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    qarray = sdpvar(prod(tempdims(:)),1);
    q = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q{coords{:}} = qarray(idx);
    end
    
    
    %q = sdpvar(nr_det_points, inputs_cell{2:end}, outputs_cell{2:end});
    
    visibility_constraints = [alpha >= 0, alpha <= 1];
    %visibility_constraints = [visibility_constraints,alpha <= 1];
    
    positivityconstraints = [];
	auxsize=size(q);
    for i=1:prod(auxsize(:))
        positivityconstraints = [positivityconstraints, qarray(i) >= 0];
    end
    
    
    %% Non signalling constraints
    
    nonsignalling_constraintsB = [];
    % non signaling for bob:
    for lam = 1:nr_det_points
        coordstructure = [nr_outputs_per_party(2), nr_inputs_per_party(2), nr_inputs_per_party(4), nr_outputs_per_party(4)];
        all_b_y_and_w_d = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_y_and_w_d,1)
            b = all_b_y_and_w_d(slice,1);
            y = all_b_y_and_w_d(slice,2);
            w = all_b_y_and_w_d(slice,3);
            d = all_b_y_and_w_d(slice,4);
            % choose z = 1 for summ1
            summ1 = 0;
            for c = 1:nr_outputs_per_party(3)
               summ1 = summ1 + q{lam,y,1,w,b,c,d};
            end
            % the marginal for z != 1 should be equal to z=1
            for z=2:nr_inputs_per_party(3)
                summ2 = 0;
                for c = 1:nr_outputs_per_party(3)
                   summ2 = summ2 + q{lam,y,z,w,b,c,d};
                end
                nonsignalling_constraintsB = [nonsignalling_constraintsB, summ1 == summ2];
            end
        end
    end

    %non signaling for charlie:
    nonsignalling_constraintsC = [];
    for lam = 1:nr_det_points
        coordstructure = [nr_outputs_per_party(3), nr_inputs_per_party(3), nr_inputs_per_party(4), nr_outputs_per_party(4)];
        all_c_z_and_w_d = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_c_z_and_w_d,1)
            c = all_c_z_and_w_d(slice,1);
            z = all_c_z_and_w_d(slice,2);
            w = all_c_z_and_w_d(slice,3);
            d = all_c_z_and_w_d(slice,4);
            % choose y = 1
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
               summ1 = summ1 + q{lam,1,z,w,b,c,d};
            end
            % the marginal for z' != 1 should be equal to z=1
            for y=2:nr_inputs_per_party(2)
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                   summ2 = summ2 + q{lam,y,z,w,b,c,d};
                end
                nonsignalling_constraintsC = [nonsignalling_constraintsC, summ1 == summ2];
            end
        end
    end
    
    %non signaling between the instrument and bob and charlie
    % nonsignaling for instrument
    % nonsignaling D not influenced by B1 and B2 (or B and C)
    nonsignalling_constraintsBCD = [];
     for lam = 1:nr_det_points
        coordstructure = [nr_outputs_per_party(4), nr_inputs_per_party(4)];
        product = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for idx = 1:size(product,1)
            d = product(idx,1);
            w = product(idx,2);
            
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
                for c = 1:nr_outputs_per_party(3)
                    summ1 = summ1 + q{lam,1,1,w,b,c,d};
                end
            end
            
            inputstructure = [nr_inputs_per_party(2), nr_inputs_per_party(3)];
            YZinputsCartesianProduct = ind2subv(inputstructure, 1:prod(inputstructure(:)));
            for index = 2:size(YZinputsCartesianProduct,1)
                y2 = YZinputsCartesianProduct(index,1);
                z2 = YZinputsCartesianProduct(index,2);
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                    for c = 1:nr_outputs_per_party(3)
                        summ2 = summ2 + q{lam,y2,z2,w,b,c,d};
                    end
                end
                nonsignalling_constraintsBCD = [nonsignalling_constraintsBCD, ...
                                        summ1 == summ2];
            end
        end
     end
     
    % partial normalization nonsignaling (summing over bcd should not depend
    % on the inputs)   
    nonsignalling_constraintsBC = [];
    for lam = 1:nr_det_points          
        inputstructure = [nr_inputs_per_party(2), nr_inputs_per_party(3), nr_inputs_per_party(4)];
        all_y_and_z_and_w = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        y1 = all_y_and_z_and_w(slice,1);
        z1 = all_y_and_z_and_w(slice,2);
        w1 = all_y_and_z_and_w(slice,3);
        summ1 = 0;
        for b = 1:nr_outputs_per_party(2)
            for c = 1:nr_outputs_per_party(3)
                for d = 1:nr_outputs_per_party(4)
                    summ1 = summ1 + q{lam,y1,z1,w1,b,c,d};
                end
            end
        end
        
        for slice = 2:size(all_y_and_z_and_w,1)
            y2 = all_y_and_z_and_w(slice,1);
            z2 = all_y_and_z_and_w(slice,2);
            w2 = all_y_and_z_and_w(slice,3);
            summ2 = 0;
            for b = 1:nr_outputs_per_party(2)
                for c = 1:nr_outputs_per_party(3)
                    for d = 1:nr_outputs_per_party(4)
                        summ2 = summ2 + q{lam,y2,z2,w2,b,c,d};
                    end
                end
            end
            nonsignalling_constraintsBC = [nonsignalling_constraintsBC, summ1 == summ2];
        end
    end
    
    %% Probability constraints
    det_strategy = givedetstratA(nr_outputs_per_party(1),nr_inputs_per_party(1));

    noisy_prob = (1-alpha) * p1 + (alpha) * p2;
    
    probability_constraints = [];
    productstructure = [nr_inputs_per_party, nr_outputs_per_party];
    cartesianproduct_forprobconstraints = ind2subv(productstructure, 1:prod(productstructure(:)));
    for idx = 1:size(cartesianproduct_forprobconstraints,1)
        coords_cell = num2cell(cartesianproduct_forprobconstraints(idx,:));
        summ = 0;
        for lam = 1:nr_det_points
            % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
            summ = summ + det_strategy(lam, coords_cell{1}, coords_cell{nrparties+1}) ... 
                * q{lam, coords_cell{2:nrparties}, coords_cell{nrparties+2:end}};
        end
        %cleanprob = clean(noisy_prob(coords_cell{:}), 1e-12);
        probability_constraints = [probability_constraints, summ == noisy_prob(coords_cell{:})];
    end
    
    %% Solving the SDP
    
    objective = alpha;
            
    constraints = [positivityconstraints, ...
        visibility_constraints, ...
        probability_constraints, ...
        nonsignalling_constraintsB, ...
        nonsignalling_constraintsC, ...
        nonsignalling_constraintsBC, ...
        nonsignalling_constraintsBCD];
    optsol = optimize(constraints, objective, ...
        sdpsettings('solver','mosek','verbose',0,'showprogress',0, 'debug', 0, 'warning',0));
    
%                       'mosek.MSK_DPAR_INTPNT_TOL_INFEAS', 1e-6,...
%                       'mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS',1e-6));
%     optsol = optimize(constraints, objective, ...
%         sdpsettings('solver','gurobi','verbose',1,'showprogress',1, 'debug', 1, 'warning',1));    
%   
    LPstatus = optsol.problem;
    % if LPstatus ~= 0
    %     disp(optsol)
    %     error('Check why the problem is not successfully solved.');
    % end
    
    dimcell = num2cell([nr_inputs_per_party,nr_outputs_per_party]);
    bellcoeffs = zeros(dimcell{:});           
    index = 1;
    for idx = 1:size(cartesianproduct_forprobconstraints,1)
        % Here I need to make sure that say, index = 6, corresponds to the
        % same (x,y,z,a,b,c) tuple that it did when defining the
        % probability constraints. I make sure of that by using the same
        % array for looping, 'cartesianproduct_forprobconstraints'
        coords_cell = num2cell(cartesianproduct_forprobconstraints(idx,:));
        bellcoeffs(coords_cell{:}) = -value(dual(probability_constraints(index)));
        index = index + 1;
    end
    %bellcoeffs = clean(bellcoeffs, 1e-16);
    final_alpha = value(alpha);
    if final_alpha>0.7
        LPstatus = 1;
        disp(final_alpha);
       error("alpha value is not good!");
    end
    dual_alpha = 0;
end
