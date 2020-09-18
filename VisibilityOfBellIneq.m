function out = VisibilityOfBellIneq(bellcoeffs, bellLocalbound, meas, channel, ins1, outs1)
if isa(ins1,'cell') || isa(outs1,'cell')
   ins = [length(ins1{1}),length(ins1{2}),length(ins1{3})];
   outs = [length(outs1{1}),length(outs1{2}),length(outs1{3})];
else
    ins = ins1;
    outs = outs1;
end


alpha = sdpvar(1,1);
alphaconstraints = [alpha >= 0, alpha <= 1];
finalstate = final_state(ini_state(alpha),channel);


cartproductOUT = ind2subv(outs, 1:prod(outs,'all'));
cartproductIN  = ind2subv(ins,  1:prod(ins,'all'));

objective = 0;
for i1 = 1:length(cartproductOUT)
    a = cartproductOUT(i1,1);
    b = cartproductOUT(i1,2);
    c = cartproductOUT(i1,3);
    for i2 = 1:length(cartproductIN)
        x = cartproductIN(i2,1);
        y = cartproductIN(i2,2);
        z = cartproductIN(i2,3);
        probability = prob(finalstate,meas,[a,b,c],[x,y,z]);
        objective = objective + bellcoeffs(x,y,z,a,b,c)*probability;
    end
end
objective = objective - bellLocalbound;

constraints = [objective <=0, alphaconstraints:'alpha'];

optsol = optimize(constraints, -objective, ...
        sdpsettings('solver','mosek', ...
    'verbose',0,'dualize',0, ...
    'showprogress',0,...
    'savesolverinput',0,'savesolveroutput',0,'debug',0));

out=value(alpha);

end

