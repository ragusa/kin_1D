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
    X_beg = X;
    
    [X,w,t,y] =  solve_prke_ode(X_beg,dt,time_end,Phi_old,Phi_new);
    p(it+1) = X(1);
end