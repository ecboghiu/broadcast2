function corr = PropToCorrelator3Parties(prob,ins,outs)
nrX = length(ins{1});
nrY = length(ins{2});
nrZ = length(ins{3});
corr = zeros(nrX + nrY + nrZ + nrX*nrY + nrX*nrZ + nrY*nrZ + nrX*nrY*nrZ);
% A_0
corr(1) = OneBodyCorrelator(prob, 1, [0,0,0], outs);
end

function oneBody = OneBodyCorrelator(prob, party, setting, outs)
    x=setting{1};
    y=setting{2};
    z=setting{3};

    if party==1
        summ0 = 0;
        summ1 = 0;
        for b=1:outs{2}
            for c=1:outs{3}
                summ0 = summ0 + prob(x,y,z,1,b,c);
                summ1 = summ1 + prob(x,y,z,2,b,c);
            end
        end
        oneBody = summ0 - summ1;
        return;
    elseif party==2
        summ0 = 0;
        summ1 = 0;
        for a=1:outs{1}
            for c=1:outs{3}
                summ0 = summ0 + prob(x,y,z,a,1,c);
                summ1 = summ1 + prob(x,y,z,a,2,c);
            end
        end
        oneBody = summ0 - summ1;
        return;
    elseif party==3
        summ0 = 0;
        summ1 = 0;
        for a=1:outs{1}
            for b=1:outs{2}
                summ0 = summ0 + prob(x,y,z,a,b,1);
                summ1 = summ1 + prob(x,y,z,a,b,2);
            end
        end
        oneBody = summ0 - summ1;
        return;
    end

end