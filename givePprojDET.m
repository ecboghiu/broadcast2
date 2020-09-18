function projs = givePprojDET()
    sig1 = [0+0i,1+0i;1+0i,0+0i];
    sig2 = [0+0i,-1i;1i,0+0i];
    sig3 = [1+0i,0+0i;0+0i,-1+0i] ;

    e1 = [1;0];
    e2 = [0;1];

    e1e1 = Tensor(e1,e1);
    e1e2 = Tensor(e1,e2);
    e2e1 = Tensor(e2,e1);
    e2e2 = Tensor(e2,e2);

    %Phi_minus = 1/2^0.5 * (e1e1 - e2e2);
    %Psi_plus  = 1/2^0.5 * (e1e2 + e2e1);
    %Phi_plus  = 1/2^0.5 * (e1e1 + e2e2);

    A0 = sig3;
    A1 = sig1;
    A2 = sig2;

    phi = atan(1/2^0.5);
    B0 = cos(phi) * sig1 + sin(phi) * sig2;
    B1 = cos(phi) * sig1 - sin(phi) * sig2;

    C0 = sig3;
    C1 = sig1;

    projs = {{giveprojs(A0),giveprojs(A1),giveprojs(A2)}, ...
                {giveprojs(B0),giveprojs(B1)}, ...
                {giveprojs(C0),giveprojs(C1)}};
end