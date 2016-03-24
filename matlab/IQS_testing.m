function [p] = IQS_testing(dt,u)

global dat npar

ntimes = length(u(1,:))-1;
time_end=0;
p = zeros(ntimes+1,1); 

[~,beff_MGT]=compute_prke_parameters(0.,u(1:npar.n,1));
X=[1;beff_MGT/dat.lambda];
p(1) = X(1);

for it=1:ntimes
    time_beg = time_end;
    time_end = time_beg+dt;
    Phi_old = u(1:npar.n,it);
    Phi_new = u(1:npar.n,it+1);
    Phi_old = u(1:npar.n,1);
    Phi_new = u(1:npar.n,1);
    X_beg = X;
    
    [X,w,t,y] =  solve_prke_ode(X_beg,dt,time_end,Phi_old,Phi_new);
    p(it+1) = X(1);
end



% save_param = npar.solve_prke_compute_rho_each_time;
% 
% npar.solve_prke_compute_rho_each_time = true;
% X2=[1;beff_MGT/dat.lambda];
% [X2,~,~,~] =  solve_prke_ode(X2,time_end,time_end,u(1:npar.n,1),u(1:npar.n,1));
% 
% error_prke_ode = (X(1)-X2(1))/X2(1)
% 
% npar.solve_prke_compute_rho_each_time = save_param;
% 
