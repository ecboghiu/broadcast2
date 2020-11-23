function corr = GiveCorrelationVector(prob, ins, outs)

corr = zeros(3+2+2+3*2+3*2+2*2+3*2*2);

% Now Ax
y=1;
z=1; % we don't care about these value for Ax because of nosignalling
for x=ins{1}
    for idx=outs{1}
        b = outs{2};
        c = outs{3};
        corr(x) = corr(x) + prob(x,y,z,idx,b,c);
    end
end

% Now By
y=1;
z=1; % we don't care about these value for Ax because of nosignalling
for x=ins{1}
    for idx=outs{1}
        b = outs{2};
        c = outs{3};
        corr(x) = corr(x) + prob(x,y,z,idx,b,c);
    end
end

% Now Cz
y=1;
z=1; % we don't care about these value for Ax because of nosignalling
for x=ins{1}
    for idx=outs{1}
        b = outs{2};
        c = outs{3};
        corr(x) = corr(x) + prob(x,y,z,idx,b,c);
    end
end

end

