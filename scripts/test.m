alpha=0.7;
state_small = ini_state(alpha);
reshaped = reshape(state_small, 2,2,2,2);

state1 = 0;
for i=1:4
    for j=1:4
        state1 =  state1 + state_small(i,j) * ketbra(i,j,4);
    end
end

state2 = 0;
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                state2 =  state2 + reshaped(i,j,k,l) * kron(ketbra(i,k,2),ketbra(j,l,2));
            end
        end    
    end
end

disp(state1-state2)
