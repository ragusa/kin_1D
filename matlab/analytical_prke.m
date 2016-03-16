function [p] = analytical_prke(u,dt,ntimes)

global dat npar

time_start  = dat.rod_mov.t_beg_1;
time_finish = dat.rod_mov.t_end_1;

it_start = time_start/dt+1;
it_end   = int16(time_finish/dt+1);
t = (0:ntimes)*dt;
p = zeros(length(t),1); rho = p;
[rho_start,~]=compute_prke_parameters(dat.rod_mov.t_beg_1,u(1:npar.n,it_start));
[rho_end  ,~]=compute_prke_parameters(dat.rod_mov.t_end_1,u(1:npar.n,it_end));

m = (rho_end-rho_start)/(time_finish-time_start);

for it=1:ntimes+1
    [rho(it),~]=compute_prke_parameters(t(it),u(1:npar.n,it));
    if it<it_start
%     Null transient part
        p(it) = 1;
    elseif it<=it_end
%     Ramp part
        p(it) = exp(m/2*(t(it)-time_start)^2);
    else
%     Constant part
        p(it) = p(it_end)*exp(rho_end*(t(it)-time_finish));
    end
end
check = t(it_end)
figure(2)
plot(t,rho)


