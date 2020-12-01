function state = ApplyMapToBipartition(inistate,channel,placement)
    % TODO put dimA, dimB, dimB1, dimB2 as function inputss
    dimstate = size(inistate,1);
    assert(mod(dimstate,2)==0,'Currently not supported to add bipartitions of uneven dimension. Should be easy to do.');
    
    dimA=dimstate/2;
    dimB=dimstate/2;
    
    dimB2=2;
    choi = ChoiMatrix(channel);
    dimB1 = size(choi,1)/dimA/2;
    dimB2 = size(choi,1)/dimA/2;

%     biggerstate = kron( inistate.', eye(dimA*dimB1*dimB2) );
%     Phi = auxPHI(dimA);
%     biggerchannel = kron(Phi*Phi',ChoiMatrix(channel));
%     swapop=Tensor(eye(dimA),SwapOperator(2),eye(dimB1*dimB2));
%     biggerchannel = swapop * biggerchannel * swapop';
%     state = PartialTrace(biggerchannel*biggerstate, [1,2], [dimA,dimB,dimA,dimB1,dimB2]);
%     

%     PartialTrace(biggerchannel, [3,4,5], [dimA,dimB,dimA,dimB1,dimB2])
%     PartialTrace(bigchoi*biggerstate, [1,2], [dimA,dimB,dimA,dimB1,dimB2])
%     PartialTrace(bigchoi, [3,4,5], [dimA,dimB,dimA,dimB1,dimB2])
%     bool1 = bigchoi == biggerchannel;
%     fprintf("is chanel ok? %d", prod(bool1(:)));
%     out1 = PartialTrace(biggerstate*biggerchannel, [1,2], [dimA,dimB,dimA,dimB1,dimB2]);
%     out2 = ApplyMap(inistate, biggerchannel);
   
    reshaped_state = reshape(inistate,dimA,dimB,dimA,dimB);
    if placement=="left"
        state = 0;
        for i=1:dimA
            for j=1:dimB
                for k=1:dimA
                    for l=1:dimB
                        state = state + reshaped_state(i,j,k,l) * ...
                                        kron( ApplyMap(ketbra(i,k,dimA),choi), ketbra(j,l,dimB));
                    end 
                end
            end
        end    
    elseif placement=="right"
        state = 0;
        for i=1:dimA
            for j=1:dimB
                for k=1:dimA
                    for l=1:dimB
                        state = state + reshaped_state(i,j,k,l) * ...
                                        kron( ketbra(i,k,dimA), ApplyMap(ketbra(j,l,dimB),choi));
                    end 
                end
            end
        end  
    else
       error("Wrong 'placement' value."); 
    end
end


