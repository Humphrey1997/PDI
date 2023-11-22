CBL1 =  importdata('CBL1.txt');
CBL2 =  importdata('CBL2.txt');
CBL3 =  importdata('CBL3.txt');
BBL1 =  importdata('BBL1.txt');
BBL2 =  importdata('BBL2.txt');
BBL3 =  importdata('BBL3.txt');

ts = tight_subplot(1,2,[.06 .08],[.14 .08],[.07 .05]);
ytext = '$\gamma /{\omega _0} $';
xtext = '$\omega /{\omega _0} $';
LHtext = '${\omega _{LH}} /{\omega _0} $';
fsize = 30;
lw = 2;

axes(ts(1));
plot(CBL1(:,2)*1e-3,CBL1(:,3)*1e-4,'Color',[0.9  0.0  0.9],'linewidth',lw)
hold on
plot(CBL2(:,2)*1e-3,CBL2(:,3)*1e-4,'Color',[0.2  0.2  1.0],'linewidth',lw)
hold on
plot(CBL3(:,2)*1e-3,CBL3(:,3)*1e-4,'Color',[0.8 0.8 0.0],'linewidth',lw)
hold on
line([0.0069 0.0069],[0 5e-4],'color','k','linestyle','--','linewidth',1.5)
hold on
text(0.0072,4.2e-4,LHtext,'color','k','FontSize',fsize,'Interpreter','latex')
text(0.0055,4.6e-4,'(a)','color','k','FontSize',fsize+10)


legend({'{\it{E}}_0=1MeV/m','{\it{E}}_0=2MeV/m','{\it{E}}_0=3MeV/m'},'Fontname','Times New Roman')
ylabel(ytext,'FontSize',fsize,'Interpreter','latex')
ylim([0 5]*1e-4)
xlim([5 15]*1e-3)
xlabel(xtext,'FontSize',fsize,'Interpreter','latex')
set(gca,'FontSize',fsize)

axes(ts(2));
plot(BBL1(:,2)*1e-3,BBL1(:,3)*1e-4,'Color',[0.9  0.0  0.9],'linewidth',lw)
hold on
plot(BBL2(:,2)*1e-3,BBL2(:,3)*1e-4,'Color',[0.2  0.2  1.0],'linewidth',lw)
hold on
plot(BBL3(:,2)*1e-3,BBL3(:,3)*1e-4,'Color',[0.8 0.8 0.0],'linewidth',lw)
hold on
line([0.0069 0.0069],[0 5e-4],'color','k','linestyle','--','linewidth',1.5)
hold on
text(0.0072,4.2e-4,LHtext,'color','k','FontSize',fsize,'Interpreter','latex')
text(0.0055,4.6e-4,'(b)','color','k','FontSize',fsize+10)

legend({'{\it{E}}_0=1MeV/m','{\it{E}}_0=2MeV/m','{\it{E}}_0=3MeV/m'},'Fontname','Times New Roman')
ylabel(ytext,'FontSize',fsize,'Interpreter','latex')
ylim([0 5]*1e-4)
xlim([5 15]*1e-3)
xlabel(xtext,'FontSize',fsize,'Interpreter','latex')
set(gca,'FontSize',fsize)

% set (gcf,'unit','normalized','position', [0,0,1,1]);
% saveas(gcf, 'test.eps');