function projs = givePprojRANDgeneral(partynr, ins, outs)
    sig1 = [0+0i,1+0i;1+0i,0+0i];
    sig2 = [0+0i,-1i;1i,0+0i];
    sig3 = [1+0i,0+0i;0+0i,-1+0i] ;
    %sig1 = Pauli(1);
    %sig2 = Pauli(2);
    %sig3 = Pauli(3);
    
    e1 = [1;0];
    e2 = [0;1];

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

%     nA0 = RandSphereSurface(3);
%     nA1 = RandSphereSurface(3);
%     nA2 = RandSphereSurface(3);
%     nA3 = RandSphereSurface(3);
%     A0 = nA0(1) * sig1 + nA0(2) * sig2 + nA0(3) * sig3;
%     A1 = nA1(1) * sig1 + nA1(2) * sig2 + nA1(3) * sig3;
%     A2 = nA2(1) * sig1 + nA2(2) * sig2 + nA2(3) * sig3;
%     A3 = nA3(1) * sig1 + nA3(2) * sig2 + nA3(3) * sig3;
% 
%     nB0 = RandSphereSurface(3);
%     nB1 = RandSphereSurface(3);
%     nB2 = RandSphereSurface(3);
%     nB3 = RandSphereSurface(3);
%     B0 = nB0(1) * sig1 + nB0(2) * sig2 + nB0(3) * sig3;
%     B1 = nB1(1) * sig1 + nB1(2) * sig2 + nB1(3) * sig3;
%     B2 = nB2(1) * sig1 + nB2(2) * sig2 + nB2(3) * sig3;
%     B3 = nB3(1) * sig1 + nB3(2) * sig2 + nB3(3) * sig3;
% 
%     nC0 = RandSphereSurface(3);
%     nC1 = RandSphereSurface(3);
%     C0 = nC0(1) * sig1 + nC0(2) * sig2 + nC0(3) * sig3;
%     C1 = nC1(1) * sig1 + nC1(2) * sig2 + nC1(3) * sig3;
%     
%     projs = {{giveprojs(A0),giveprojs(A1),giveprojs(A2),giveprojs(A3)}, ...
%                 {giveprojs(B0),giveprojs(B1),giveprojs(B2),giveprojs(B3)}, ...
%                 {giveprojs(C0),giveprojs(C1)}};
end