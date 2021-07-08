function change_of_basis_mat = fromBellToCorrMat(ins, outs)
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
% We should ignore the element (1,1,1) as that has no correspondence in 
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
cell_dims = num2cell(dims);

corrdims = ins + 1; % This is because we take into accout the identity as well, (Id, A1, A2) etc.
cell_corrdims = num2cell(corrdims);

% Matrix where each column corresponds to a correlator and the column is
% the Bell coefficients of the correlator, namely its expression in terms
% of p(abc|xyz)
change_of_basis_mat = zeros(prod(corrdims), prod(dims)); 


inputs_slices = ind2subv(ins, 1:prod(ins(:)));
outputs_slices = ind2subv(outs, 1:prod(outs(:)));
for idx_ins = 1:size(inputs_slices,1) % Ignore (1,1,1) as that is a constant, the local bound, so up to size(corr_coords,1)-1
    coords_in = num2cell(inputs_slices(idx_ins, :));
    
    x = coords_in{1};
    y = coords_in{2};
    z = coords_in{3};
    
    for idx_outs = 1:size(outputs_slices,1)
        coords_out = num2cell(outputs_slices(idx_outs, :));

        a = (-1)^(coords_out{1}-1);
        b = (-1)^(coords_out{2}-1);
        c = (-1)^(coords_out{3}-1);
        
        correlator = zeros(cell_corrdims{:});
        correlator(1,1,1) = 1;
        correlator(x+1,1,1) = a;
        correlator(1,y+1,1) = b;
        correlator(1,1,z+1) = c;
        correlator(x+1,y+1,1) = a*b;
        correlator(x+1,1,z+1) = a*c;
        correlator(1,y+1,z+1) = b*c;
        correlator(x+1,y+1,z+1) = a*b*c;
        correlator = correlator / 8;
        
        indice = sub2ind([ins, outs], coords_in{:}, coords_out{:});
        
        change_of_basis_mat(:, indice) = correlator(:);
    end
end

end