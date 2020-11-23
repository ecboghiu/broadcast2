function U = give_SingletForBC()
    e1 = [1;0];
    e2 = [0;1];

    e1e1 = Tensor(e1,e1);
    e1e2 = Tensor(e1,e2);
    e2e1 = Tensor(e2,e1);
    e2e2 = Tensor(e2,e2);

    Phi_minus = 1/2^0.5 * (e1e1 - e2e2);
    Psi_plus  = 1/2^0.5 * (e1e2 + e2e1);

    alp  = pi/8;
    psi0 =  sin(alp) * Phi_minus + cos(alp) * Psi_plus;
    psi1 = -cos(alp) * Phi_minus + sin(alp) * Psi_plus;
    
    U = Phi_minus * e1';
end
