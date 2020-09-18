function state = ini_state(vis)
    e1 = [1;0];
    e2 = [0;1];

    e1e1 = Tensor(e1,e1);
    %e1e2 = Tensor(e1,e2);
    %e2e1 = Tensor(e2,e1);
    e2e2 = Tensor(e2,e2);

    % Phi_minus = 1/2^0.5 * (e1e1 - e2e2);
    % Psi_plus  = 1/2^0.5 * (e1e2 + e2e1);
    Phi_plus  = 1/2^0.5 * (e1e1 + e2e2);
    
    state = vis * (Phi_plus * Phi_plus') + (1-vis) * eye(4)/4;
end

