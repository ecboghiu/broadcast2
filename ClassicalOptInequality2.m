function funcoutput = ClassicalOptInequality2(bellcoeffs,ins1, outs1) 
    % ins is something like {[1,2,3],[1,2],[1,2]} giving the 
    % inputs of each party. cell 1 is for party A, cell 2 is for party B
    % and cell 3 is for party C. Similarly outs is something of the
    % form {[1,2],[1,2],[1,2,3,4]} each cell element giving the different
    % outputs of each party.
    % in practice we only care about the NUMBER of inputs or outputs
    % and never about what they are exactly
    
    % TODO Fix this
    ins = {1:ins1(1),1:ins1(2),1:ins1(3)};
    outs = {1:outs1(1),1:outs1(2),1:outs1(3)};
    
    party_for_det_points = 1; % this is party 'A'
    nrinputsofA = length(ins{party_for_det_points});
    nroutputsofA = length(outs{party_for_det_points});
    nr_det_points = 1;
    for i=1:nrinputsofA
       nr_det_points = nr_det_points*nroutputsofA; 
    end
 
%     alpha = sdpvar(1); % alpha will be the visibility
%     finalstate = final_state(ini_state(alpha),channel);
%     visibility_constraints = [alpha >= 0, alpha <= 1];
    
    q = sdpvar(nr_det_points*length(ins{2})*length(ins{3})* ...
                             length(outs{2})*length(outs{3}),1);
                         

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
       for y = ins{2}
           for z = ins{3}
               for b = outs{2}
                   for c = outs{3}
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
        [Out,In] = meshgrid(outs{2},ins{2});
        prod = [Out(:) In(:)];
        for idx = 1:length(prod(:,1))
            b = prod(idx,1);
            y = prod(idx,2);
            combinations = nchoosek(ins{3},2);
            for idx2 = 1:length(combinations(:,1))
                z1 = combinations(idx2,1);
                z2 = combinations(idx2,2);
                summ1 = 0;
                for idx3 = outs{3}
                    summ1 = summ1 + qmatrix{lam}{y}{z1}{b}{idx3};
                end
                
                summ2 = 0;
                for idx3 = outs{3}
                    summ2 = summ2 + qmatrix{lam}{y}{z2}{b}{idx3};
                end
                nonsignalling_constraints = [nonsignalling_constraints, ...
                                            summ1 == summ2];
            end
        end
    end
    %non signaling for charlie:
    for lam = 1:nr_det_points
        [Out,In] = meshgrid(outs{3},ins{3});
        prod = [Out(:) In(:)];
        for idx = 1:length(prod(:,1))
            c = prod(idx,1);
            z = prod(idx,2);

            combinations = nchoosek(ins{2},2);
            for idx2 = 1:length(combinations(:,1))
                y1 = combinations(idx2,1);
                y2 = combinations(idx2,2);
                summ1 = 0;
                for idx3 = outs{2}
                    summ1 = summ1 + qmatrix{lam}{y1}{z}{idx3}{c};
                end
                summ2 = 0;
                for idx3 = outs{2}
                    summ2 = summ2 + qmatrix{lam}{y2}{z}{idx3}{c};
                end
                nonsignalling_constraints = [nonsignalling_constraints, ...
                                        summ1 == summ2];
            end
        end
    end

    % probability constraints
    det = givedetstratA(length(outs{1}),length(ins{1}));


%     probability_constraints = [];
%     probability_constraints_inp_out = [];
%     
%     [X,Y,Z] = meshgrid(outs{1},outs{2},outs{3});
%     cartproductOUT = [X(:) Y(:) Z(:)];
%     [X,Y,Z] = meshgrid(ins{1},ins{2},ins{3});
%     cartproductIN = [X(:) Y(:) Z(:)];
%     for i1 = 1:length(cartproductOUT)
%        for i2 = 1:length(cartproductIN)
%           a = cartproductOUT(i1,1);
%           b = cartproductOUT(i1,2);
%           c = cartproductOUT(i1,3);
%           x = cartproductIN(i2,1);
%           y = cartproductIN(i2,2);
%           z = cartproductIN(i2,3);
% 
%           summ = 0;
%           for aux = 1:nr_det_points
%               summ = summ + det(aux,x,a) * qmatrix{aux}{y}{z}{b}{c};
%           end
% 
%           probability = prob(finalstate,meas,[a,b,c],[x,y,z]);
%           probability_constraints = [probability_constraints, ...
%                                         summ == probability];
%           probability_constraints_inp_out = ...
%               [probability_constraints_inp_out, [a;b;c;x;y;z]];
%        end
%     end

    
    [X,Y,Z] = meshgrid(outs{1},outs{2},outs{3});
    cartproductOUT = [X(:) Y(:) Z(:)];
    [X,Y,Z] = meshgrid(ins{1},ins{2},ins{3});
    cartproductIN = [X(:) Y(:) Z(:)];
    
    objective = 0;
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
          objective = objective + bellcoeffs(x,y,z,a,b,c)*summ;
       end
    end
    
    constraints = [...%probability_constraints:'probability', ...
                ...%visibility_constraints:'visibility', ...
                positivityconstraints(1):'pos1', ...
                positivityconstraints(2:end):'positivity', ...
                nonsignalling_constraints:'nonsignalling'];

    optsol = optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
    'verbose',0,'dualize',0, ...
    'showprogress',0,...
    'savesolverinput',0,'savesolveroutput',0,'debug',0));

% uncomment if you have probablity constraints and you want the bell
% inequality
%     nrduals = length(probability_constraints);
%     dualvals = zeros(nrduals,1);
%     for i=1:nrduals
%         dualvals(i) = value(dual(probability_constraints(i)));
%     end

                
    funcoutput = value(objective);
end

