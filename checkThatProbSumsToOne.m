function checkThatProbSumsToOne(prob,ins,outs)

tol = 1e-5;
    for x=1:ins(1)
        for y=1:ins(2)
            for z=1:ins(3)
                summ=0;
                for a=1:outs(1)
                    for b=1:outs(2)
                        for c=1:outs(3)
                        summ = summ + prob(x,y,z,a,b,c);
                    end
                end
            end
            if abs(summ -1)>tol
                fprintf("Something went wrong. sum_outputs(prob) = %f, (x,y,z)=(%d,%d,%d)", summ, x, y, z);
            end
        end
    end
end
                        
                        

end
