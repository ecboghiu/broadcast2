function summ = dispBellCoeffsINSTR(bellcoeffs,ins,outs)
    summ = 0;
    structure = [ins,outs];
    loopvars = ind2subv(structure, 1:prod(structure,'all'));
    for idx = 1:size(loopvars,1)
        ascell = num2cell(loopvars(idx,:));
        x = loopvars(idx,1);
        y = loopvars(idx,2);
        z = loopvars(idx,3);
        w = loopvars(idx,4);
        a = loopvars(idx,5);
        b = loopvars(idx,6);
        c = loopvars(idx,7);
        d = loopvars(idx,8);
        str = join(['p_',string(a),...
                         string(b),...
                         string(c),...
                         string(d),...
                         '_',...
                         string(x),...
                         string(y),...
                         string(z),...
                         string(w)],'');
        var = sym(char(str));
        summ = summ + var*bellcoeffs(ascell{:});
    end
    summ = vpa(summ,3);
end

