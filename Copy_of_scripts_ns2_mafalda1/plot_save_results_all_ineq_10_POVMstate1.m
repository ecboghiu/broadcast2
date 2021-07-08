load('results_all_ineq_10_POVMstate.mat');

CONST_CHI = 0.02;
COS2 = (cos(2*CONST_CHI))^2;
roots23 = roots([COS2,-2*COS2,0,2,-1]);
roots23 = roots23(abs(imag(roots23))<1e-8); % discard imaginary roots;
roots23 = roots23(abs(roots23)>=0);
roots23 = roots23(abs(roots23)<=1); % discard p outside [0,1]
roots23 = max(roots23); % just in case there is more than one, but there shouldnt be
p = roots23 - 0.1; % check that ineq (23) is satisfied in the intervat (0,roots23)
assert(dot([COS2,-2*COS2,0,2,-1],[p^4,p^3,p^2,p^1,1]) <= 0, "Something bad happened");
% for all p >
UNSTEERABILITY_THRESHOLD = roots23; % we're using opposite convention for noise

figure
plot(1-[results_per_ineq{:,1}])
xlim([0, NR_OF_INEQS+1])
%ylim([0.8,1])

set(gca,'XTick',0:4:NR_OF_INEQS);

yline(UNSTEERABILITY_THRESHOLD,'-','','LineWidth',1,'DisplayName','p_threshold');
