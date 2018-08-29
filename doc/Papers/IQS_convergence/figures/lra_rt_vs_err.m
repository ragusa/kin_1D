clear; clc; close all

a = 0.5;

bf.err = [1.407e-2 3.174e-3 7.690e-4 1.892e-4 4.590e-5];
bf.rt = [4.11 6.01 10.38 21.91 25.23];
bf.dt = [4e-3, 2e-3, 1e-3, 0.5e-3, 0.25e-3];
bf.eff = bf.err.^a.*bf.rt;

iqs.err = [2.612e-3 9.893e-4 5.796e-4 4.772e-4 4.516e-4];
iqs.rt = [3.96 6.02 7.87 12.61 22.14];
iqs.nt = [1, 2, 4, 8, 16];
iqs.eff = iqs.err.^a.*iqs.rt;

iqspc.err = [3.488e-3 1.349e-3 9.161e-4 8.052e-4 7.905e-4];
iqspc.rt = [2.91 3.73 3.97 5.39 8.19];
iqspc.nt = [1, 2, 4, 8, 16];
iqspc.eff = iqspc.err.^a.*iqspc.rt;

figure
loglog(iqs.err,iqs.rt,'o-', iqspc.err,iqspc.rt,'o-', bf.err,bf.rt,'o-')
xlabel('Error')
ylabel('Runtime (hr)')
legend('IQS, dt=0.004, N_T=1, 2, 4, 8, 16','IQS P-C, dt=0.004, N_T=1, 2, 4, 8, 16','Implicit Dis., dt=4e-3, 2e-3, 1e-3, 0.5e-3, 0.25e-3')
grid on

figure
semilogy(iqs.nt,iqs.eff,'o-', iqspc.nt,iqspc.eff,'o-')
hold on
nn = linspace(min(iqs.nt),max(iqs.nt),100);
for i=1:length(bf.eff)
    eff = ones(size(nn))*bf.eff(i);
    semilogy(nn,eff,'--');
end

