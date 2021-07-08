function newPOVMs= eta_POVMs(povms, party, strat, eta, ins, outs)

assert(0 <= eta <= 1, "Invalid value for eta.");

newPOVMs = povms;

ze = zeros(outs(party));
id = eye(outs(party));

if strat==1
    det = {{id,ze},{id,ze}}; % x=0->a=0, x=1->a=0
elseif strat==2
    det = {{id,ze},{ze,id}}; % x=0->a=0, x=1->a=1
elseif strat==3
    det = {{ze,id},{id,ze}}; % x=0->a=1, x=1->a=0
elseif strat==4
    det = {{ze,id},{ze,id}}; % x=0->a=1, x=1->a=1
else
    error("invalid strategy");
end
    
for x = 1:ins(party)
   for a=1:outs(party)
       newPOVMs{party}{x}{a} = eta * povms{party}{x}{a} + (1-eta)* det{x}{a};
   end
end


end

