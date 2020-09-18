function summ = dispBellCoeffs(bellcoeffs,ins,outs)
    summ = 0;
    for x=ins{1}
        for y=ins{2}
            for z=ins{3}
                for a=outs{1}
                    for b=outs{2}
                        for c=outs{3}
                            str = join(['p_',string(a),...
                                             string(b),...
                                             string(c),...
                                             '_',...
                                             string(x),...
                                             string(y),...
                                             string(z)],'');
                            var = sym(char(str));
                            summ = summ + var*bellcoeffs(x,y,z,a,b,c);
                        end
                    end
                end
            end
        end
    end
    summ = vpa(summ,3);
end

