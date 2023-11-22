EM3B =  importdata('EM3B.txt');
EM3C =  importdata('EM3C.txt');
EM4B =  importdata('EM4B.txt');
EM4C =  importdata('EM4C.txt');

ts = tight_subplot(1,2,[.06 .08],[.14 .08],[.07 .05]);
ytext = '$\gamma /{\omega _0} $';
xtext = '$\omega /{\omega _0} $';
LHtext = '${\omega _{LH}} /{\omega _0} $';
fsize = 30;
lw = 2;

axes(ts(1));
plot(EM4C(:,2)*1e-3,EM4C(:,3)*1e-4,'Color',[0.9  0.0  0.9],'linewidth',lw)
hold on
plot(EM4B(:,2)*1e-3,EM4B(:,3)*1e-4,'Color',[0.2  0.2  1.0],'linestyle','--','linewidth',lw)
hold on
text(0.0072,4.6e-4,'(a)','color','k','FontSize',fsize+10)

legend({'X-mode pump','EBW-L pump'},'Fontname','Times New Roman')
ylabel(ytext,'FontSize',fsize,'Interpreter','latex')
ylim([0 5]*1e-4)
xlim([7 11]*1e-3)
set(gca,'XTick',0.007:0.002:0.011);
set(gca,'XTickLabel',{'0.007','0.009','0.011'})
xlabel(xtext,'FontSize',fsize,'Interpreter','latex')
set(gca,'FontSize',fsize)

axes(ts(2));
plot(EM3C(:,2)*1e-3,EM3C(:,3)*1e-4,'Color',[0.9  0.0  0.9],'linewidth',lw)
hold on
plot(EM3B(:,2)*1e-3,EM3B(:,3)*1e-4,'Color',[0.2  0.2  1.0],'linestyle','--','linewidth',lw)
hold on
text(0.0072,4.6e-4,'(b)','color','k','FontSize',fsize+10)

legend({'X-mode pump','EBW-L pump'},'Fontname','Times New Roman')
ylabel(ytext,'FontSize',fsize,'Interpreter','latex')
ylim([0 5]*1e-4)
xlim([7 11]*1e-3)
set(gca,'XTick',0.007:0.002:0.011);
set(gca,'XTickLabel',{'0.007','0.009','0.011'})
xlabel(xtext,'FontSize',fsize,'Interpreter','latex')
set(gca,'FontSize',fsize)
