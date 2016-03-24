function [p,p_matlab] = analytical_prke(u,dt,ntimes)

global dat npar

time_start  = dat.rod_mov.t_beg_1;
time_finish = dat.rod_mov.t_end_1;

it_start = time_start/dt+1;
it_end   = int16(time_finish/dt+1);
t = (0:ntimes)*dt;
p = zeros(length(t),1); rho = p;
[rho_start,~]=compute_prke_parameters(0   ,u(1:npar.n,1));
[rho_end  ,~]=compute_prke_parameters(1e20,u(1:npar.n,1));

m = (rho_end-rho_start)/(time_finish-time_start);

for it=1:ntimes+1
%     [rho(it),~]=compute_prke_parameters(t(it),u(1:npar.n,it));
    if it<=it_start
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


h = @(t,t1,t2) 0*(t<=t1)+m*(t-t1).*(t<=t2).*(t>t1)+m*(t2-t1)*(t>t2);

% tolerances for odesolvers
rtol = 3e-14; atol = 3e-14;
options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);
% call ode solver with prke function that re-compute rho/beta at each time (expensive)
[t_matlab,p_matlab]=ode15s(@(t,y) h(t,time_start,time_finish)*y,[0 ntimes*dt],1.,options);

check2 = (p(end)-p_matlab(end))/p_matlab(end)
