data_FS1 = importdata('zfs1.txt');
data_FS2 = importdata('zfs2.txt');
data_FS3 = importdata('zfs3.txt');

%不同kz的最大增长率
lw = 2;
fsize = 30;
wci = 2.4216*1e-4;
veci = 208;

ts = tight_subplot(1,2,[.06 .09],[.14 .08],[.07 .03]);
axes(ts(1));
plot(data_FS1(:,1)*veci,data_FS1(:,3)*1e-4/wci,'--r','linewidth',lw)
hold on
plot(data_FS2(:,1)*veci,data_FS2(:,3)*1e-4/wci,':b','linewidth',lw)
hold on
plot(data_FS3(:,1)*veci,data_FS3(:,3)*1e-4/wci,'-.m','linewidth',lw)
hold on
text(0.4,5.5,'(a)','color','k','FontSize',fsize+10)

legend('N=1','N=2','N=3','location','NorthEast','FontSize',fsize)
xlim([0 5])
ylim([0 6])
xtext = '${k_{sz}}{v_e}/{\omega _{ci}}$';
ytext = '${\omega}/{\omega _{ci}}$';
xlabel(xtext,'FontSize',fsize,'Interpreter','latex')
ylabel(ytext,'FontSize',fsize,'Interpreter','latex')
set(gca,'fontsize',fsize)


axes(ts(2));
plot(data_FS1(:,1)*veci,data_FS1(:,4)*1e-6/wci,'--r','linewidth',lw)
hold on
plot(data_FS2(:,1)*veci,data_FS2(:,4)*1e-6/wci,':b','linewidth',lw)
hold on
plot(data_FS3(:,1)*veci,data_FS3(:,4)*1e-6/wci,'-.m','linewidth',lw)
hold on
text(0.4,5.5/6*0.024,'(b)','color','k','FontSize',fsize+10)

legend('N=1','N=2','N=3','location','NorthEast','FontSize',fsize)
xlim([0 5])
ylim([0 0.024])
xtext = '${k_{sz}}{v_e}/{\omega _{ci}}$';
ytext = '$\gamma /{\omega _{ci}} $';
xlabel(xtext,'FontSize',fsize,'Interpreter','latex')
ylabel(ytext,'FontSize',fsize,'Interpreter','latex')
ax = gca;
ax.YAxis.Exponent = -2;
set(gca,'fontsize',fsize)

