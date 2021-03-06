function [final_alpha,bellcoeffs,LPstatus,dual_alpha] = BroadcastNS2_LP(p1,p2)
    % IMPORTANT alpha is defined as (1-alpha)*p1 + (alpha)*p2
    % usually p1 is the entangled, and p2 the uniform, so this
    % program minimizes alpha using this convention
    % if p1 and p2 are both local then all alpha are good
    % if p1 and p2 are both nonlocal the program is infeasible
    % ideally p1 should be nonlocal

    % IMPORTANT I use the convention where the probabilities are called as
    % p(x,y,z,a,b,c) instead of p(a,b,c,x,y,z)
    dims_p = size(p1);
    
    assert(mod(length(dims_p),2) == 0, "The probability array should have equal number of inputs and outputs.");
    assert(all(size(p1) == size(p2)), "The two probability arrays should have equal dimenisons.");
    assert( norm(p1(:) - p2(:)) > 1e-6, "The probabilities should not be almost equal."); 
    
    nrparties = length(dims_p)/2;
    
    nr_inputs_per_party = dims_p(1:nrparties);
    nr_outputs_per_party = dims_p(nrparties+1:end);


    %nr_det_points = nroutputsofA^nrinputsofA;
    nr_det_points = nr_outputs_per_party.^nr_inputs_per_party;
  
    alpha = sdpvar(1); % alpha will be the visibility

    %% For lam
    tempdims = [nr_det_points(3), nr_inputs_per_party(1:2), nr_outputs_per_party(1:2)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    lamarray = sdpvar(prod(tempdims(:)),1);
    q_lam = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q_lam{coords{:}} = lamarray(idx);
    end
    
    %% For mu
    tempdims = [nr_det_points(2), nr_inputs_per_party(1), nr_inputs_per_party(3), nr_outputs_per_party(1), nr_outputs_per_party(3)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    muarray = sdpvar(prod(tempdims(:)),1);
    q_mu = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q_mu{coords{:}} = muarray(idx);
    end
    
    %% For nu
    tempdims = [nr_det_points(1), nr_inputs_per_party(2:end), nr_outputs_per_party(2:end)];
    q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
    tempdims_cell = num2cell(tempdims);
    nuarray = sdpvar(prod(tempdims(:)),1);
    q_nu = cell(tempdims_cell{:});
    for idx = 1:size(q_coords,1)
        coords = num2cell(q_coords(idx,:));
        q_nu{coords{:}} = nuarray(idx);
    end
    

    visibility_constraints = [alpha >= 0];

    positivityconstraints = [];
	auxsize=size(q_lam);
    for i=1:prod(auxsize(:))
        positivityconstraints = [positivityconstraints, lamarray(i) >= 0];
    end
    auxsize=size(q_mu);
    for i=1:prod(auxsize(:))
        positivityconstraints = [positivityconstraints, muarray(i) >= 0];
    end
    auxsize=size(q_nu);
    for i=1:prod(auxsize(:))
        positivityconstraints = [positivityconstraints, nuarray(i) >= 0];
    end
    
    
    %% Non signalling constraints
    %% For q_nu
    nonsignalling_constraintsB_nu = [];
    % non signaling for bob:
    for nu = 1:nr_det_points(1)
        coordstructure = [nr_outputs_per_party(2), nr_inputs_per_party(2)];
        all_b_and_y = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_and_y,1)
            b = all_b_and_y(slice,1);
            y = all_b_and_y(slice,2);
            % choose z = 1
            summ1 = 0;
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q_nu{nu,y,1,b,c};
            end
            % the marginal for z != 1 should be equal to z=1
            for z=2:nr_inputs_per_party(3)
                summ2 = 0;
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q_nu{nu,y,z,b,c};
                end
                nonsignalling_constraintsB_nu = [nonsignalling_constraintsB_nu, summ1 == summ2];
            end
        end
    end
    %non signaling for charlie:
    nonsignalling_constraintsC_nu = [];
    for nu = 1:nr_det_points(1)
        coordstructure = [nr_outputs_per_party(3), nr_inputs_per_party(3)];
        all_b_and_y = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_and_y,1)
            c = all_b_and_y(slice,1);
            z = all_b_and_y(slice,2);
            % choose y = 1
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
                summ1 = summ1 + q_nu{nu,1,z,b,c};
            end
            % the marginal for z' != 1 should be equal to z=1
            for y=2:nr_inputs_per_party(2)
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                    summ2 = summ2 + q_nu{nu,y,z,b,c};
                end
                nonsignalling_constraintsC_nu = [nonsignalling_constraintsC_nu, summ1 == summ2];
            end
        end
    end
    % partial normalization nonsignaling (summing over bc should not depend
    % on the inputs)   
    nonsignalling_constraintsBC = [];
    for nu = 1:nr_det_points(1)          
        inputstructure = [nr_inputs_per_party(2), nr_inputs_per_party(3)];
        all_y_and_z = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        y1 = all_y_and_z(slice,1);
        z1 = all_y_and_z(slice,2);
        summ1 = 0;
        for b = 1:nr_outputs_per_party(2)
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q_nu{nu,y1,z1,b,c};
            end
        end
        for slice = 2:size(all_y_and_z,1)
            y2 = all_y_and_z(slice,1);
            z2 = all_y_and_z(slice,2);
            summ2 = 0;
            for b = 1:nr_outputs_per_party(2)
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q_nu{nu,y2,z2,b,c};
                end
            end
            nonsignalling_constraintsBC = [nonsignalling_constraintsBC, summ1 == summ2];
        end
    end
    
    %% For mu
    nonsignalling_constraintsA_mu = [];
    % non signaling for alice:
    for mu = 1:nr_det_points(2)
        coordstructure = [nr_outputs_per_party(1), nr_inputs_per_party(1)];
        all_a_and_x = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_a_and_x,1)
            a = all_a_and_x(slice,1);
            x = all_a_and_x(slice,2);
            % choose z = 1
            summ1 = 0;
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q_mu{mu,x,1,a,c};
            end
            % the marginal for z != 1 should be equal to z=1
            for z=2:nr_inputs_per_party(3)
                summ2 = 0;
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q_mu{mu,x,z,a,c};
                end
                nonsignalling_constraintsA_mu = [nonsignalling_constraintsA_mu, summ1 == summ2];
            end
        end
    end
    %non signaling for charlie:
    nonsignalling_constraintsC_mu = [];
    for mu = 1:nr_det_points(2)
        coordstructure = [nr_outputs_per_party(3), nr_inputs_per_party(3)];
        all_c_and_z = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_c_and_z,1)
            c = all_c_and_z(slice,1);
            z = all_c_and_z(slice,2);
            % choose x = 1
            summ1 = 0;
            for a = 1:nr_outputs_per_party(1)
                summ1 = summ1 + q_mu{mu,1,z,a,c};
            end
            % the marginal for z' != 1 should be equal to z=1
            for x=2:nr_inputs_per_party(1)
                summ2 = 0;
                for a = 1:nr_outputs_per_party(1)
                    summ2 = summ2 + q_mu{mu,x,z,a,c};
                end
                nonsignalling_constraintsC_mu = [nonsignalling_constraintsC_mu, summ1 == summ2];
            end
        end
    end
    % partial normalization nonsignaling (summing over bc should not depend
    % on the inputs)   
    nonsignalling_constraintsAC = [];
    for mu = 1:nr_det_points(2)          
        inputstructure = [nr_inputs_per_party(1), nr_inputs_per_party(3)];
        all_x_and_z = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        x1 = all_x_and_z(slice,1);
        z1 = all_x_and_z(slice,2);
        summ1 = 0;
        for a = 1:nr_outputs_per_party(1)
            for c = 1:nr_outputs_per_party(3)
                summ1 = summ1 + q_mu{mu,x1,z1,a,c};
            end
        end
        
        for slice = 2:size(all_x_and_z,1)
            x2 = all_x_and_z(slice,1);
            z2 = all_x_and_z(slice,2);
            summ2 = 0;
            for a = 1:nr_outputs_per_party(1)
                for c = 1:nr_outputs_per_party(3)
                    summ2 = summ2 + q_mu{mu,x2,z2,a,c};
                end
            end
            nonsignalling_constraintsAC = [nonsignalling_constraintsAC, summ1 == summ2];
        end
    end
    
    %% For lam
    nonsignalling_constraintsA_lam = [];
    % non signaling for alice:
    for lam = 1:nr_det_points(3)
        coordstructure = [nr_outputs_per_party(1), nr_inputs_per_party(1)];
        all_a_and_x = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_a_and_x,1)
            a = all_a_and_x(slice,1);
            x = all_a_and_x(slice,2);
            % choose y = 1
            summ1 = 0;
            for b = 1:nr_outputs_per_party(2)
                summ1 = summ1 + q_lam{lam,x,1,a,b};
            end
            % the marginal for y != 1 should be equal to y=1
            for y=2:nr_inputs_per_party(2)
                summ2 = 0;
                for b = 1:nr_outputs_per_party(2)
                    summ2 = summ2 + q_lam{lam,x,y,a,b};
                end
                nonsignalling_constraintsA_lam = [nonsignalling_constraintsA_lam, summ1 == summ2];
            end
        end
    end
    %non signaling for bob:
    nonsignalling_constraintsB_lam = [];
    for lam = 1:nr_det_points(3)
        coordstructure = [nr_outputs_per_party(2), nr_inputs_per_party(2)];
        all_b_and_y = ind2subv(coordstructure, 1:prod(coordstructure(:)));
        for slice = 1:size(all_b_and_y,1)
            b = all_b_and_y(slice,1);
            y = all_b_and_y(slice,2);
            % choose y = 1
            summ1 = 0;
            for a = 1:nr_outputs_per_party(1)
                summ1 = summ1 + q_lam{lam,1,y,a,b};
            end
            % the marginal for z' != 1 should be equal to z=1
            for x=2:nr_inputs_per_party(1)
                summ2 = 0;
                for a = 1:nr_outputs_per_party(1)
                    summ2 = summ2 + q_lam{lam,x,y,a,b};
                end
                nonsignalling_constraintsB_lam = [nonsignalling_constraintsB_lam, summ1 == summ2];
            end
        end
    end
    % partial normalization nonsignaling (summing over bc should not depend
    % on the inputs)   
    nonsignalling_constraintsAB = [];
    for lam = 1:nr_det_points(3)          
        inputstructure = [nr_inputs_per_party(1), nr_inputs_per_party(2)];
        all_x_and_y = ind2subv(inputstructure, 1:prod(inputstructure(:)));
        
        slice = 1;
        x1 = all_x_and_y(slice,1);
        z1 = all_x_and_y(slice,2);
        summ1 = 0;
        for a = 1:nr_outputs_per_party(1)
            for b = 1:nr_outputs_per_party(2)
                summ1 = summ1 + q_lam{lam,x1,z1,a,b};
            end
        end
        
        for slice = 2:size(all_x_and_y,1)
            x2 = all_x_and_y(slice,1);
            z2 = all_x_and_y(slice,2);
            summ2 = 0;
            for a = 1:nr_outputs_per_party(1)
                for b = 1:nr_outputs_per_party(2)
                    summ2 = summ2 + q_lam{lam,x2,z2,a,b};
                end
            end
            nonsignalling_constraintsAB = [nonsignalling_constraintsAB, summ1 == summ2];
        end
    end
    
    
    %% Probability constraints
    det_strategyA = givedetstratA(nr_outputs_per_party(1),nr_inputs_per_party(1));
    det_strategyB = givedetstratA(nr_outputs_per_party(2),nr_inputs_per_party(2));
    det_strategyC = givedetstratA(nr_outputs_per_party(3),nr_inputs_per_party(3));

    noisy_prob = (1-alpha) * p1 + (alpha) * p2;
    
    probability_constraints = [];
    productstructure = [nr_inputs_per_party, nr_outputs_per_party];
    cartesianproduct_forprobconstraints = ind2subv(productstructure, 1:prod(productstructure(:)));
    for idx = 1:size(cartesianproduct_forprobconstraints,1)
        coords_cell = num2cell(cartesianproduct_forprobconstraints(idx,:));
        summ = 0;
        for nu = 1:nr_det_points(1)
            % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
            summ = summ + det_strategyA(nu, coords_cell{1}, coords_cell{nrparties+1}) ... 
                * q_nu{nu, coords_cell{2:nrparties}, coords_cell{nrparties+2:end}};
        end
        for mu = 1:nr_det_points(2)
            % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
            summ = summ + det_strategyB(mu, coords_cell{2}, coords_cell{nrparties+2}) ... 
                * q_mu{mu, coords_cell{1}, coords_cell{3}, coords_cell{nrparties+1}, coords_cell{nrparties+3}};
        end
        for lam = 1:nr_det_points(3)
            % for 3 parties, same as det_strategy(lam,x,a) * q(lam,y,z,b,c);
            summ = summ + det_strategyC(lam, coords_cell{3}, coords_cell{nrparties+3}) ... 
                * q_lam{lam, coords_cell{1:2}, coords_cell{nrparties+1:nrparties+2}};
        end
        
        %cleanprob = clean(noisy_prob(coords_cell{:}), 1e-12);
        probability_constraints = [probability_constraints, summ == noisy_prob(coords_cell{:})];
    end
    
    %% Solving the SDP
    
    objective = alpha;
            
    constraints = [positivityconstraints, ...
        visibility_constraints, ...
        probability_constraints, ...
        nonsignalling_constraintsA_lam, ...
        nonsignalling_constraintsA_mu, ...
        nonsignalling_constraintsB_lam, ...
        nonsignalling_constraintsB_nu, ...
        nonsignalling_constraintsC_nu, ...
        nonsignalling_constraintsC_mu, ...
        nonsignalling_constraintsBC, ...
        nonsignalling_constraintsAC, ...
        nonsignalling_constraintsAB];
    
    optsol = optimize(constraints, objective, ...
        sdpsettings('solver','mosek','verbose',0,'showprogress',0, 'debug', 0, 'warning',0));
    
%                       'mosek.MSK_DPAR_INTPNT_TOL_INFEAS', 1e-6,...
%                       'mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS',1e-6));
%     optsol = optimize(constraints, objective, ...
%         sdpsettings('solver','gurobi','verbose',1,'showprogress',1, 'debug', 1, 'warning',1));    
%   
    LPstatus = optsol.problem;
%     if LPstatus ~= 0
%         disp(optsol)
%         error('Check why the problem is not successfully solved.');
%     end
    
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
    final_alpha = value(alpha);
    dual_alpha = 0;
    %fprintf("lam4 · (p1-p2) = %f\n", sum(bellcoeffs .* (p1-p2),'all'));
end
