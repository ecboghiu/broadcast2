function change_of_basis_mat = fromCorrToBellMat_INSTR(ins, outs)
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
probs_shape = zeros(celldims{:});  % !! coeffs called as probs_shape(x,y,z,w,a,b,c,d) instead of probs_shape(a,b,c,d,x,y,z,w)

% !!! this is because in the instrumental scenario we have signalling from the instrument
% so we have things like <AxByw> etc.
% Notice that in <Ax Byw Czw Dw> actually all the w's are the same
% To encode this correlator I use (x,y,z,w,w2) with x,y,z,w taking values
% up to the ones specified in "ins" but also the value 0 (1 in matlab)
% which encodes the Identity. The value w2 is a trigger showing wether the
% correlator Dw appears or not. While for x, x=0, denotes Ax=Id, we cannot
% use w=0 do denote Dw=Id because then we cannot asign a value do Byw and
% Czw. And as such we have w always take a value in (1, max(w)=ins(4)) and
% have a separate flag for when Dw appears.
% numerically, it will be a matrix of the sizes specified as follows.
corrdims = [1+ins(1:3), 2, ins(4)]; 

% Matrix where each column corresponds to a correlator and the column is
% the Bell coefficients of the correlator, namely its expression in terms
% of p(abc|xyz)
change_of_basis_mat = zeros(prod(corrdims), numel(probs_shape)); 

outputs_slices = ind2subv(outs, 1:prod(outs));
inputs_slices = ind2subv(ins, 1:prod(ins));
corr_coords = ind2subv(corrdims, 1:prod(corrdims(:)));
nr_coords = numel(corr_coords(1,:));

corrdims_cell = num2cell(corrdims);
aux_correlator_array = zeros(corrdims_cell{:});
prob_size_idx = 1;
for out_idx = 1:size(outputs_slices,1)
    outputs = outputs_slices(out_idx,:);
    %outputs_slice = num2cell(outputs);
    for in_idx = 1:size(inputs_slices,1)
       inputs_slice = num2cell(inputs_slices(in_idx,:));
       
       aux_correlator_array = 0*aux_correlator_array;
       for corr_idx = 1:(size(corr_coords,1))  % first element(1,1,1) corresponds to a constant
                                               % TODO What about w=2?
           coords = corr_coords(corr_idx,:);
           coordscell = num2cell(coords); % this containts w as well, coords = [x y z w2 w] so actually 2 entries correspond to constant value
           x = coords(1);
           y = coords(2);
           z = coords(3);
           w2 = coords(4);
           w = coords(5);
 
           if all([x,y,z,w2,w] == [1, 1, 1, 1, 1])
               aux_correlator_array(coordscell{:}) = 1; % multiply by 1/16 at the end
           else
               % find which of [x,y,z,w2] are different from 1
               % note that w is never different from 1 in this convention
               % remember that everything is displaced by 1, so if x,y,z
               % are 1 this means Ax=Id, By=Id, Cz=Id
               x = coords(1);
               y = coords(2);
               z = coords(3);
               w2 = coords(4);
               w = coords(5);
               
               which_are_not_1 = [x,y,z,w2]~=1 ;
               % Now from the ones which are not 1 compare with (x,y,z,w)
               % form prob(abcd|xyzw)
               ins_corr = [[x,y,z]-1,w];
               ins_prob = [inputs_slice{:}];
               if ins_prob(which_are_not_1) ~= ins_corr(which_are_not_1)
                  aux_correlator_array(coordscell{:}) = 0;
               else
                  % if not, then set the coefficient
                  % example: abd<AxBywDw>
                  % we need to set aux_correlator_array(x,y,1,w2=1,w)=abd.
                  aux_correlator_array(coordscell{:}) = prod((-1).^(outputs(which_are_not_1)-1));
               end       
           end
       end
       change_of_basis_mat(:, prob_size_idx) = 1/16*aux_correlator_array(:);
       aux_correlator_array = 0*aux_correlator_array;
       
       prob_size_idx = prob_size_idx + 1;
    end
end

end