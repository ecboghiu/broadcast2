function flag=checkThatProbSumsToOne(prob,ins,outs)
flag = true;
tol  = 1e-3;

cartproductOUT = ind2subv(outs, 1:prod(outs,'all'));
cartproductIN  = ind2subv(ins,  1:prod(ins,'all'));
for i2 = 1:length(cartproductIN)
    inputs = cartproductIN(i2,:);
    summ = 0;
    for i1 = 1:length(cartproductOUT)
        outputs = cartproductOUT(i1,:);
        indexes = num2cell([inputs,outputs]);
        summ = summ + prob(indexes{:});
    end
    if abs(summ-1)>tol
        fprintf("Something went wrong. sum_outputs(prob) = %f, (x,y,z)=",summ);
        disp(inputs);
        flag = false;
        return
    end
end        
                        
end

