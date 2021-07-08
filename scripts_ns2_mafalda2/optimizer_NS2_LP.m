function opt = optimizer_NS2_LP(ins, outs)
    nrparties = length(ins);
    
    nr_inputs_per_party = ins;
    nr_outputs_per_party = outs;
    
    auxdims = num2cell([ins,outs]);
    prob1 = sdpvar(auxdims{:},'full','real');
    prob2 = sdpvar(auxdims{:},'full','real'); 


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

    coordstructure = [ins, outs];
    auxdims = num2cell([ins,outs]);
    noisy_prob = cell(auxdims{:});
    all_a_and_x = ind2subv(coordstructure, 1:prod(coordstructure(:)));
    for slice = 1:size(all_a_and_x,1)
        auxdims = num2cell(all_a_and_x(slice, :));
        noisy_prob{auxdims{:}} = (1-alpha) * prob1(auxdims{:}) + (alpha) * prob2(auxdims{:});
    end
    %noisy_prob = (1-alpha) * prob1 + (alpha) * prob2;
    
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
        
        probability_constraints = [probability_constraints, summ == noisy_prob{coords_cell{:}}];
    end
    
    %% Solving the SDP
    
    objective = alpha;
            
    constraints = [probability_constraints, ...
        positivityconstraints, ...
        visibility_constraints, ...
        nonsignalling_constraintsA_lam, ...
        nonsignalling_constraintsA_mu, ...
        nonsignalling_constraintsB_lam, ...
        nonsignalling_constraintsB_nu, ...
        nonsignalling_constraintsC_nu, ...
        nonsignalling_constraintsC_mu, ...
        nonsignalling_constraintsBC, ...
        nonsignalling_constraintsAC, ...
        nonsignalling_constraintsAB];
    
    opt = cell(2);
    nr_duals = length(probability_constraints);
    opt{1} = optimizer(constraints, objective, ...
        sdpsettings('solver','mosek'), {prob1, prob2}, objective);
    opt{2} = nr_duals; % for extracting the LP
    
end