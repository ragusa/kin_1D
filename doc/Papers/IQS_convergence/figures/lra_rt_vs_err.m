clear; clc; close all

bf.err = [1.407e-2 3.174e-3 7.690e-4 1.892e-4 4.590e-5];
bf.rt = [4.11 6.01 10.38 21.91 25.23];
bf.dt = [4e-3, 2e-3, 1e-3, 0.5e-3, 0.25e-3];
bf.eff = bf.err.*bf.rt;

iqs.err = [2.612e-3 9.893e-4 5.796e-4 4.772e-4 4.516e-4];
iqs.rt = [3.96 6.02 7.87 12.61 22.14];
iqs.nt = [1, 2, 4, 8, 16];
iqs.eff = iqs.err.*iqs.rt;

iqspc.err = [3.488e-3 1.349e-3 9.161e-4 8.052e-4 7.905e-4];
iqspc.rt = [2.91 3.73 3.97 5.39 8.19];
iqspc.nt = [1, 2, 4, 8, 16];
iqspc.eff = iqspc.err.*iqspc.rt;

figure
lw = 1.5;
ms = 6;
loglog(iqs.err,iqs.rt,'o-','LineWidth',lw,'MarkerSize',ms)
hold on
loglog(iqspc.err,iqspc.rt,'o-','LineWidth',lw,'MarkerSize',ms)
loglog(bf.err,bf.rt,'o-','LineWidth',lw,'MarkerSize',ms)
set(gca,'YLim',[2.91 40])
xlabel('Error','Interpreter','latex','FontSize',14)
ylabel('Run-time (hr)','Interpreter','latex','FontSize',14)
leg = {'IQS, $\Delta t=$0.004, $N_T=$1, 2, 4, 8, 16',...
       'IQS-PC, $\Delta t=$0.004, $N_T=$1, 2, 4, 8, 16',...
       'Implicit Dis., $\Delta t$=4e-3, 2e-3, 1e-3, 0.5e-3, 0.25e-3'};
legend(leg(:),'Interpreter','latex','FontSize',12)
grid on

figure
loglog(0.004./iqs.nt,iqs.eff,'o-',0.004./iqspc.nt,iqspc.eff,'o-',bf.dt,bf.eff,'o-')
set(gca,'XTick',bf.dt(end:-1:1))
grid on
% line(iqs.nt,iqs.eff,'Color','b','Marker','o')
% hold on
% line(iqspc.nt,iqspc.eff,'Color','b','Marker','^')
% ax1 = gca;
% ax1.XColor = 'b';
% ax1.YScale = 'log';
% ax1.YLim = [1e-3 0.1];
% ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YScale','log','Color','none');
% line(bf.dt,bf.eff,'Color','k','Marker','s','Parent',ax2)
% hold on
% nn = linspace(0,max(iqs.nt),100);
% for i=1:length(bf.eff)
%     eff = ones(size(nn))*bf.eff(i);
%     semilogy(nn,eff,'--');
% end
% grid on

