% THIS DOESNT WORK AS OF NOW
function out = auxSwap(input, dims, chosen2)
% important: the output is not the swap operator, but the swap operator
% applied to the input
swapdim_1 = chosen2(1);
swapdim_2 = chosen2(2);
if length(chosen2)>2 || length(chosen2) <2
   error('can only swap two subsystems'); 
end

celldims = num2cell([dims,dims]);
reshapedInput = reshape(input, celldims{:});

newdims = dims;
temp = newdims(swapdim_1);
newdims(swapdim_1) = dims(swapdim_2);
newdims(swapdim_2) = temp;
newcelldims = num2cell([newdims,newdims]);
out = zeros(newcelldims{:}); 

coordstructure = [dims, dims];
cartesianproduct = ind2subv(coordstructure, 1:prod(coordstructure(:)));
for idx = 1:size(cartesianproduct,1)
    coords = num2cell(cartesianproduct(idx,:));
    coords_new = coords;
  
    % working thorugh an example, say coords is {i,j,k,l,m,n}
    % and swapdim_1 = 2, swampdim_2 = 3
    % then coords_new we want it to be
    % {i, k, j, l, n, m}
    % see it in braket form
    % M = \sum M_ijklmn |ijk><lmn| = \sum M_ijklmn |i><l| \otimes |j><l| \otimes |k><n|
    % after swapping |j><l| with |k><n| we get the desired output
    temp = coords_new{swapdim_1};
    coords_new{swapdim_1} = coords{swapdim_2};
    coords_new{swapdim_2} = temp;
    offset = length(dims);
    temp = coords_new{offset+swapdim_1};
    coords_new{offset+swapdim_1} = coords{offset+swapdim_2};
    coords_new{offset+swapdim_2} = temp;
    %disp([coords,coords_new])
    
    % now map the values of the operator onto the new one
    %disp(reshapedInput(coords{:}))
    out(coords_new{:}) = reshapedInput(coords{:});
end
    
out = reshape(out, size(input)); % make it a square matrix again
end
