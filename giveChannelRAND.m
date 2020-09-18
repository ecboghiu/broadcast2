function channel = giveChannelRAND(dimx,dimy)
    id = eye(dimx);
    
    e = cell(dimx);
    for i=1:dimx
       e{i} = id(:,i); 
    end
    
    o = cell(dimy);
    channel = 0;
    for i=1:dimx
       o{i} = (-1+2*rand(dimy,1))+1i*(-1+2*rand(dimy,1));
       o{i} = o{i}/norm(o{i});
       channel = channel + o{i} * e{i}.';
    end
end

