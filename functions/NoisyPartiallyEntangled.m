function state = NoisyPartiallyEntangled(p,xi)
    e1 = [1;0];
    e2 = [0;1];

    e1e1 = Tensor(e1,e1);
    e1e2 = Tensor(e1,e2);
    e2e1 = Tensor(e2,e1);
    e2e2 = Tensor(e2,e2);
    
    psi_xi = cos(xi) * e1e1 + sin(xi) * e2e2;
    psi_xi = psi_xi * psi_xi'/ norm(psi_xi);
    
    rhoA = PartialTrace(psi_xi, 2, [2,2]);
    rhoB = PartialTrace(psi_xi, 1, [2,2]);
    
    state = (1-p) * psi_xi + p * kron(eye(2)/2, rhoB);
    %state = (1-p) * psi_xi + p * kron(rhoA, eye(2)/2);
end

