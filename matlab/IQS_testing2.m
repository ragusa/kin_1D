function [p] = IQS_testing2(dt,ntimes,u)

global dat npar io

time_end=0;
p = zeros(ntimes+1,1);

[~,beff_MGT]=compute_prke_parameters(0.,u(1:npar.n,1));
X=[1;beff_MGT/dat.lambda];
p(1) = X(1);

npar.iqs_prke_interpolation_method = 4;

for it=1:ntimes

    time_beg = time_end;
    time_end = time_beg+dt;
    dat.ode.time_beg=time_beg;
    dat.ode.time_end=time_end;

    Phi_old = io.phi_save(:,it);
    Phi_new = io.phi_save(:,it+1);
    ode_beg = io.odefsave(:,it);
    ode_end = io.odefsave(:,it+1);

    X_beg = X;
    
    % tolerances for odesolvers
    rtol = 3e-14; abso = 3e-14;
    atol  = abso*ones(length(X),1);
    options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);
    
    t1=time_beg;
    t2=time_end;
    mat=[ t1^3 t1^2 t1 1; t2^3 t2^2 t2 1;  3*t1^2 2*t1 1 0; 3*t2^2 2*t2 1 0];
    rhs = [ Phi_old' ;Phi_new' ; ode_beg'; ode_end'];
    % contains the w coefficients such that :
    %  shape(t) = w(1) t^3 + w(2) t^2 + w(3) t + w(4)
    w = mat\rhs; w=w'; 
    dat.ode.nw = size(w,2);
    for k=1:dat.ode.nw
        [dat.ode.rho_MGT_beg(k),dat.ode.beff_MGT_beg(k)]=compute_prke_parameters(time_beg,w(:,k));
        [dat.ode.rho_MGT_end(k),dat.ode.beff_MGT_end(k)]=compute_prke_parameters(time_end,w(:,k));
    end
    
    % call ode solver with prke function that linearly interpolates
    [t,y]=ode15s(@funprke_linhermite_interp,[time_beg time_end],X,options);
    
    % [X,w,t,y] =  solve_prke_ode(X_beg,dt,time_end,Phi_old,Phi_new);
    X=y(:,ned);
    p(it+1) = X(1);

end


