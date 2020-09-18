function checkThatProbSumsToOne(prob,ins,outs)

tol = 1e-5;

for x=ins{1}
    for y=ins{2}
        for z=ins{3}
            summ = 0;
            for a=outs{1}
                for b=outs{2}
                    for c=outs{3}
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

