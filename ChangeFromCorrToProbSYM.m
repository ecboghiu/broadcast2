function prob_expression = ChangeFromCorrToProbSYM(corr_expression, ins, outs)

% THE UGLY WAY

for x=ins{1}
    for y=ins{2}
        for z=ins{3}
            for a=outs{1}
                for b=outs{2}
                    for c=outs{3}
                        var = ToCorrelatorNotation([x-1,y-1,z-1],[a-1,b-1,c-1]);
                        summ = summ + var*bellcoeffs(x,y,z,a,b,c);
                    end
                end
            end
        end
    end
end

end

