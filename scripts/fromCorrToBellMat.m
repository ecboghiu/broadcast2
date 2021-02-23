function change_of_basis_mat = fromCorrToBellMat(ins, outs)
% We want a change of basis matrix from a vector of numbers in the basis
% [A1, A2, B1, A1*B1, A2*B1, B2, A1*B2, A2*B2, ...
%                   C1, A1*C1, A2*C1, B1*C1, A1*B1*C1, A1*B1*C1, B2*C1, A1*B2*C1, A2*B2*C1, ...
%                   C2, A1*C2, A2*C2, B1*C2, A1*B1*C2, A2*B1*C2, B2*C2, A1*B2*C2, A2*B2*C2];
% To something in the basis [p(abc|xyz)]. Note that this change of basis
% is not unique.
% The pattern is to look at the product of (Id, A1, A2) x (Id, B1, B2) x
% (Id, C1, C2). MATLAB starts counting from left to right meaning it is
% 000 100 010 110 
% instead of the traditional 
% 000, 001, 010, 011, etc.
% This happens when using ind2subv.
% We should ignore the element (1,1,1) as that has no correspondente in 
% p(abc|xyz).

% For the MATLAB behaviour of counting from left to right instead of right
% to left uncomment the following to convince yourself
% % % tempdims = [ins, outs];
% % % dims = num2cell(tempdims);
% % % aux = 1:prod(tempdims(:));
% % % aux_mat = reshape(aux, dims{:});  % the reshape uses the built-in matlab counting
% % % q_coords = ind2subv(tempdims, 1:prod(tempdims(:)));
% % % for idx = 1:size(q_coords,1)
% % %     coords = num2cell(q_coords(idx,:));
% % %     assert(aux_mat(coords{:})==aux(idx));

nr_parties = length(ins);

dims = [ins, outs];
celldims = num2cell(dims);
% preassign memory for faster computation (?)
bellcoeffs = zeros(celldims{:});  % !! coeffs called as bellcoeffs(x,y,z,a,b,c) instead of bellcoeffs(a,b,c,x,y,z)


corrdims = ins + 1; % This is because we take into accout the identity as well, (Id, A1, A2) etc.

% Matrix where each column corresponds to a correlator and the column is
% the Bell coefficients of the correlator, namely its expression in terms
% of p(abc|xyz)
change_of_basis_mat = zeros( numel(bellcoeffs), prod(corrdims)-1); 

corr_coords = ind2subv(corrdims, 1:prod(corrdims(:)));
nr_coords = numel(corr_coords(1,:));
for idx = 1:(size(corr_coords,1)-1)  % Ignore (1,1,1) as that is a constant, the local bound, so up to size(corr_coords,1)-1
    coords = corr_coords(idx+1,:);
    coordscell = num2cell(coords);
    %if ~all(coords == ones(1,nr_coords)) 
        outputs_slices = ind2subv(outs, 1:prod(outs));
        for out_idx = 1:size(outputs_slices,1) 
            outputs_slice = num2cell(outputs_slices(out_idx,:));
            inputs_slice = coords - 1;
            bell_coefficient = 1;
            for aux_in = 1:size(inputs_slice,2)
               if inputs_slice(aux_in) ~= 0
                   disp(bell_coefficient)
                  bell_coefficient = bell_coefficient * (-1)^(outputs_slice{aux_in}-1);
               end
            end
            for aux_idx = 1:size(inputs_slice,2)
               if inputs_slice(aux_idx) == 0
                   % To understand this suppose we have A1*B2. Then this
                   % corresponds to coords=(2,3,1). We substract 1
                   % to get the corresponding x=1 and y=2, but then the 
                   % elements which are 0 correspond to z=(whatever). Since
                   % we assume no-signalling, then we don't care about
                   % party C. We will just choose a random input for it,
                   % which I choose to be z=1 just for consistency (we can 
                   % also throw a random number but I feel this might give
                   % some problems if it changes every time).
                   inputs_slice(aux_idx) = 1;
               end
            end
            inputs_slice_cell = num2cell(inputs_slice);
            bellcoeffs(inputs_slice_cell{:},outputs_slice{:}) = bell_coefficient;
        end
        change_of_basis_mat(:, idx) = bellcoeffs(:);
    %end
    
    bellcoeffs = 0*bellcoeffs;
end
end