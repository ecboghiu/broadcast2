function [summ, scaling] = dispBellCoeffsCorrelators(bellcoeffs,ins,outs)
    scaling=1;
    summ = 0;
    
    cartproductOUT = ind2subv(outs, 1:prod(outs,'all'));
    cartproductIN  = ind2subv(ins,  1:prod(ins, 'all'));
    for i1 = 1:size(cartproductOUT,1)
        a = cartproductOUT(i1,1);
        b = cartproductOUT(i1,2);
        c = cartproductOUT(i1,3);
        for i2 = 1:size(cartproductIN,1)
            x = cartproductIN(i2,1);
            y = cartproductIN(i2,2);
            z = cartproductIN(i2,3);
            var = ToCorrelatorNotationWithInOut([x,y,z],[a,b,c]);
            summ = summ + var*bellcoeffs(x,y,z,a,b,c);
        end
    end
    summ = expand(simplify(summ));
    
    % clean the terms
    [C,T] = coeffs(summ);
    tolerance = 1e-6;
    C(abs(C)<tolerance)=0;
    summ = dot(C,T);
    
    % normalize by smallest one
    [C,T] = coeffs(summ); % redo the split so there are no 0 coeffs to divide by
    scaling = 1.0/min(abs(C(:)));
    C = scaling*C;
    
    summ = dot(C,T);
    
    summ = vpa(summ,3); % vpa to simplify integer fractions
end

