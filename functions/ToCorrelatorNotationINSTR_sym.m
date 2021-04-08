function sym_out = ToCorrelatorNotationINSTR_sym(correlator_vec, ins, outs)

corrdims = [1+ins(1:3), 2, ins(4)]; 

correlator = reshape(correlator_vec, corrdims);

parties = ['A','B','C','D'];

corr_coords = ind2subv(corrdims, 1:prod(corrdims(:)));
summ = 0;
for corr_idx = 1:(size(corr_coords,1))  % first element(1,1,1) corresponds to a constant
    coords = corr_coords(corr_idx,:);
    coordscell = num2cell(coords); 
    x = coords(1);
    y = coords(2);
    z = coords(3);
    w2 = coords(4);
    w = coords(5);
    inputs = [x,y,z,w2];
    which_are_not_1 = ~(inputs==1);
    
    included_parties = parties(which_are_not_1);
    inputs2 = [[x,y,z]-1, w];
    included_inputs = inputs2(which_are_not_1);
    assert(length(included_parties)==length(included_inputs),"Something wrong");
    
    producto = 1;
    for idx=1:length(included_parties)
        if included_parties(idx) == 'A'
            Ax = join([included_parties(idx),string(included_inputs(idx)),''],''); 
        	Ax = sym(char(Ax));
            producto = producto * Ax;
        elseif included_parties(idx) == 'B'
            Byw = join([included_parties(idx),string(included_inputs(idx)),string(included_inputs(end)),''],''); 
        	Byw = sym(char(Byw));
            producto = producto * Byw;
        elseif included_parties(idx) == 'C'
            Czw = join([included_parties(idx),string(included_inputs(idx)),string(included_inputs(end)),''],''); 
        	Czw = sym(char(Czw));
            producto = producto * Czw;
        elseif included_parties(idx) == 'd'
            Dw = join([included_parties(idx),string(included_inputs(end)),''],''); 
        	Dw = sym(char(Dw));
            producto = producto * Dw;
        end

    end
    summ = summ + correlator(coordscell{:})*producto;
end
sym_out = summ;
end

