function projs = givePprojRANDgeneral(partynr, ins, ~)
    sig1 = [0+0i,1+0i;1+0i,0+0i];
    sig2 = [0+0i,-1i;1i,0+0i];
    sig3 = [1+0i,0+0i;0+0i,-1+0i] ;
    %sig1 = Pauli(1);
    %sig2 = Pauli(2);
    %sig3 = Pauli(3);
    
%     e1 = [1;0];
%     e2 = [0;1];

    %e1e1 = Tensor(e1,e1);
    %e1e2 = Tensor(e1,e2);
    %e2e1 = Tensor(e2,e1);
    %e2e2 = Tensor(e2,e2);

    %Phi_minus = 1/2^0.5 * (e1e1 - e2e2);
    %Psi_plus  = 1/2^0.5 * (e1e2 + e2e1);
    %Phi_plus  = 1/2^0.5 * (e1e1 + e2e2);
    
    projs={{0}};
    for party = 1:partynr
       for pIn = 1:ins(party)
           nP = RandSphereSurface(3);
           P = nP(1) * sig1 + nP(2) * sig2 + nP(3) * sig3;
           projs{party}{pIn} = giveprojs(P); % built-in with outs=[2,2,2], only qubits
       end
    end
end