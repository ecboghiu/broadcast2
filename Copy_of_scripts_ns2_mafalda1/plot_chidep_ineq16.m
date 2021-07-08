clear all

load('chi_dependence_ineq.mat');
%load('matlab_plot_chi_ineq2.mat');

h = figure;

chiresults{end-3,2} = 0.2692801673;
chiresults{end-4,2} = 0.292829;
chiresults{end-2,2} = 0.292847;
chiresults{end-1,2} = 0.2928216117;
chiresults{2,2} = 0.292807;

chiresults{1,2} = 0.283005;

% 0.001000 0.278457 0.012519

xdata = [chiresults{:,1}];
ydata = [chiresults{:,2}];
zdata = [chiresults{:,3}];

yline(1/sqrt(2),'-','','LineWidth',1,'DisplayName','1/√2 line');


hold on
plot(xdata,1-ydata,'DisplayName','p_{critical}','LineWidth',1.8)
plot(xdata,1-zdata,'DisplayName','p_{threshold}','LineWidth',1)

title('NS2 inequality nr. 2 visibility depending on $\chi$', 'Interpreter', 'latex');
xlim([0 1.6]);
ylim([0 1]);
xlabel('$\chi$','Interpreter', 'latex')
ylabel('Visibility $p$', 'Interpreter', 'latex')

xticks([0 0.3 0.6 0.7854 0.9 1.2 1.5, 1.5708, 1.8]);
xticklabels({'','','','π/4','','','','π/2',''});
legend('Location','southeast')

annotation('textbox', [0.175 0.215 0.215 0.], 'interpreter','latex','String','State: $\rho=p \vert \psi_\chi \rangle \langle \psi_\chi \vert + (1-p) \frac{\mathrm{1}}{2}\otimes \rho^B_\chi$','FitBoxToText','on')


grid on
