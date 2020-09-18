function state = final_state(inistate,channel)
    % Note that this is not general it assumes some structure already

    dimA=2;
    dimB=2;
    
    reshaped_state = reshape(inistate,dimA,dimB,dimA,dimB);
    
    state = 0;
    for i=1:dimA
        for j=1:dimB
            for k=1:dimA
                for l=1:dimB
                    state = state + reshaped_state(i,j,k,l) * ...
              kron( ketbra(i,k,dimA), ApplyMap(ketbra(j,l,dimB),channel));
                end 
            end
        end
    end
end
