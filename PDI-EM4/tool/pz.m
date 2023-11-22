DR =  importdata('res.txt');
h = figure(1);
lcolor = [0.9 0.1 0.1];
rcolor = [0.35 0.6 0.5];
set(h,'defaultAxesColorOrder',[lcolor;rcolor])
yyaxis left
plot(DR(:,1),DR(:,3),'Color',lcolor,'linewidth',2)
hold on
% ylim([0 20])
yltext = '$\gamma /{\omega _0} {10^{ - 6}}$';
ylabel(yltext,'FontSize',20,'Interpreter','latex')
yyaxis right
plot(DR(:,1),DR(:,2),'Color',rcolor,'linewidth',2)
yrtext = '${\omega_r} /{\omega _0} {10^{ - 4}}$';
ylabel(yrtext,'FontSize',20,'Interpreter','latex')
% ylim([0 10])



xtext = '${n_ \bot } = {k_ \bot }c/{\omega _0}$';
xlabel(xtext,'FontSize',30,'Interpreter','latex')
% title('Nonlinear Dispersion Relation($e{\phi _0}/T = 0.02$,$n_0z = 0.1$,$n_Lz = 0.01$)','FontSize',20,'Interpreter','latex')
set(gca,'FontSize',30)