function [visibility, LPstatus] = visibilityOfBellInequality(bellcoeffs, localbound, p_entangled, p_uniform)

assert(mod(length(size(bellcoeffs)),2)==0,"There should be as many inputs as outputs.");
dims = size(bellcoeffs);
nrparties = length(dims)/2;
ins = dims(1:nrparties);
outs = dims(nrparties+1:end);

alpha = sdpvar(1);
vis_constraints = [alpha>=0, alpha<=1];

p_noisy = (1-alpha)*p_entangled + alpha*p_uniform;
% objective = sum(bellcoeffs .* p_noisy, size(bellcoeffs));
% aux = bellcoeffs .* p_noisy;
% auxsize=size(aux);
% while prod(auxsize(:)) ~= 1
%    aux = sum(aux); 
% end

objective = 0;
auxsize=size(bellcoeffs);
coords = ind2subv(size(bellcoeffs), 1:prod(auxsize(:)));
for idx = 1:size(coords,1)
    coords_choice = num2cell(coords(idx,:));
    objective = objective + bellcoeffs(coords_choice{:}) * p_noisy(coords_choice{:});
end
%sdisplay(objective);
constraints = [vis_constraints, objective <= localbound];

optsol = optimize(constraints, -objective, sdpsettings('solver','mosek','verbose', 0));
LPstatus = optsol.problem;
if optsol.problem ~= 0
    disp('error in visibility LP');
   %error('check why this LP is not working'); 
end

visibility = value(alpha);
end

