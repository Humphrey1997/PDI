ES3_1 =  importdata('ES3_0.4.txt');
ES3_2 =  importdata('ES3_0.8.txt');
ES3_3 =  importdata('ES3_1.2.txt');
ES4_1 =  importdata('ES4_0.4.txt');
ES4_2 =  importdata('ES4_0.8.txt');
ES4_3 =  importdata('ES4_1.2.txt');

ts = tight_subplot(1,2,[.06 .09],[.14 .08],[.09 .03]);
ytext = '$\gamma /{\omega _{ci}} $';
xtext = '$\omega /{\omega _{ci}} $';
wci = 2.4216*1e-4;
fsize = 30;
lw = 2;

axes(ts(1));
plot(ES4_1(:,2)*1e-4/wci,ES4_1(:,3)*1e-6/wci,'--r','linewidth',lw)
hold on
plot(ES4_2(:,2)*1e-4/wci,ES4_2(:,3)*1e-6/wci,':b','linewidth',lw)
hold on
plot(ES4_3(:,2)*1e-4/wci,ES4_3(:,3)*1e-6/wci,'-.m','linewidth',lw)
hold on
text(1.958,5.5/6*0.024,'(a)','color','k','FontSize',fsize+10)

legend({'\zeta = 0.4','\zeta = 0.8','\zeta = 1.2'},'Fontname','Times New Roman')
ylabel(ytext,'FontSize',fsize,'Interpreter','latex')
ylim([0 0.024])
xlim([1.95 2.05])
xlabel(xtext,'FontSize',fsize,'Interpreter','latex')
ax = gca;
ax.YAxis.Exponent = -2;
set(gca,'FontSize',fsize)

axes(ts(2));
plot(ES3_1(:,2)*1e-4/wci,ES3_1(:,3)*1e-6/wci,'--r','linewidth',lw)
hold on
plot(ES3_2(:,2)*1e-4/wci,ES3_2(:,3)*1e-6/wci,':b','linewidth',lw)
hold on
plot(ES3_3(:,2)*1e-4/wci,ES3_3(:,3)*1e-6/wci,'-.m','linewidth',lw)
hold on
text(1.958,5.5/6*0.024,'(b)','color','k','FontSize',fsize+10)

legend({'\zeta = 0.4','\zeta = 0.8','\zeta = 1.2'},'Fontname','Times New Roman')
ylabel(ytext,'FontSize',fsize,'Interpreter','latex')
ylim([0 0.024])
xlim([1.95 2.05])
xlabel(xtext,'FontSize',fsize,'Interpreter','latex')
ax = gca;
ax.YAxis.Exponent = -2;
set(gca,'FontSize',fsize)
