function p = giveProbNDArray(state, Pproj, ins, outs)
dims = num2cell([ins outs]);
p = zeros(dims{:});

cartproductOUT = ind2subv(outs, 1:prod(outs(:)));
cartproductIN  = ind2subv(ins,  1:prod(ins(:)));
for i1 = 1:length(cartproductOUT)
    a = cartproductOUT(i1,1);
    b = cartproductOUT(i1,2);
    c = cartproductOUT(i1,3);
    for i2 = 1:length(cartproductIN)
        x = cartproductIN(i2,1);
        y = cartproductIN(i2,2);
        z = cartproductIN(i2,3);
        p(x,y,z,a,b,c) = prob(state,Pproj,[a,b,c],[x,y,z]);
    end
end
end

