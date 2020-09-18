CartesianProduct({[1,2],[2],[1,2,3]})

function cart = CartesianProduct(inputcell)
%UNTITLED2 Does the cartesian product of the columns of the input
%   Detailed explanation goes here

nrcells = length(inputcell);

sizeofcells = [];
for i=1:nrcells
    sizeofcells = [sizeofcells, length(inputcell{i})];
    disp(length(inputcell{i}));
end

idx=zeros(length(sizeofcells));
for i=1:length(sizeofcells)
    for idx(i)=1:2
        disp(idx(i));
    end
    end
    
end

end

