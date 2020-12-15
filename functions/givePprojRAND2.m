function projs = givePprojRAND2(ins,outs)
nrparties = length(ins);
assert(length(ins)==length(outs),'not equal vec lengths');

projs={};
for party=1:nrparties
    dim_sys = outs(party);
    auxcell = {};
    for x=1:ins(party)
       auxcell{x} = RandomPOVM(dim_sys,outs(party));
    end
    projs{party} = auxcell;
end
%     projs = {{RandomPOVM(2,2),RandomPOVM(2,2),RandomPOVM(2,2)}, ...
%                 {RandomPOVM(2,2),RandomPOVM(2,2)}, ...
%                 {RandomPOVM(2,2),RandomPOVM(2,2)}};
end