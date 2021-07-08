
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
newdir2 = strcat(newdir,filesep,'functions',filesep);
addpath(newdir);
addpath(newdir2);

load('results_all_ineq_10_POVMstate.mat');

for ineq_nr=1:NR_OF_INEQS

COS2 = (cos(2*CONST_CHI))^2;
roots23 = roots([COS2,-2*COS2,0,2,-1]);
roots23 = roots23(abs(imag(roots23))<1e-8); % discard imaginary roots;
roots23 = roots23(abs(roots23)>=0);
roots23 = roots23(abs(roots23)<=1); % discard p outside [0,1]
roots23 = max(roots23); % just in case there is more than one, but there shouldnt be
p = roots23 - 0.1; % check that ineq (23) is satisfied in the intervat (0,roots23)
assert(dot([COS2,-2*COS2,0,2,-1],[p^4,p^3,p^2,p^1,1]) <= 0, "Something bad happened");
% for all p >
p_threshold = roots23; % we're using opposite convention for noise
   
bellcoeffs = bellcoeffs_cell{ineq_nr};
localboundNS2 = local_upper_bounds(ineq_nr);
quantumbound = table3arXiv11122626(ineq_nr, 5);

channel = results_per_ineq{ineq_nr, 2};
POVMs   = results_per_ineq{ineq_nr, 3};

p_entangled = ProbMultidimArray(final_state(PartiallyEntangledPOVM(0, CONST_CHI, 'A'), channel), POVMs, ins, outs);
p_uniform   = ProbMultidimArray(final_state(PartiallyEntangledPOVM(1, CONST_CHI, 'A'), channel), POVMs, ins, outs);

b1 = bellcoeffs.*p_entangled; b1 = sum(b1(:));
b2 = bellcoeffs.*p_uniform; b2 = sum(b2(:));

[p_crit, ~]= visibilityOfBellInequality(bellcoeffs, localboundNS2, p_entangled, p_uniform);

if p_crit == 0
fprintf("ineq nr: %d \tbell(p=0) NS2bound bell(p=1) = \t%.10f \t%f \t%.10f \tp_crit=%f \tp_thresh=%f (no violation)\n", ineq_nr, p_crit, 1-p_threshold, b1, localboundNS2, b2); 
else
if p_crit > 1-p_threshold
fprintf("ineq nr: %d \tbell(p=0) NS2bound bell(p=1) = \t%.10f \t%f \t%.10f \tp_crit=%f \tp_thresh=%f (good! p_crit > p_thresh)\n", ineq_nr, b1, localboundNS2, b2 , p_crit, 1-p_threshold); 
else
fprintf("ineq nr: %d \tbell(p=0) NS2bound bell(p=1) = \t%.10f \t%f \t%.10f \tp_crit=%f \tp_thresh=%f \n", ineq_nr, b1, localboundNS2, b2 , p_crit, 1-p_threshold); 
end
end

end