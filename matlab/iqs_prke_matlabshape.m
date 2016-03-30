function [p] = iqs_prke_matlabshape(dt,ntimes,t_ref,exact_arr)

global dat npar

time_end=0;
p = zeros(ntimes+1,1);

[~,beff_MGT]=compute_prke_parameters(0.,exact_arr(1:npar.n,1));
X=[1;beff_MGT/dat.lambda];
p(1) = X(1);

% % for it=1:ntimes
% % 
% %     time_beg = time_end;
% %     time_end = time_beg+dt;
% %     X_beg = X;
% %     
% %     % tolerances for odesolvers
% %     rtol = 3e-14; abso = 3e-14;
% %     atol  = abso*ones(length(X),1);
% %     options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);    
% %     
% %     % call ode solver with prke function that linearly interpolates
% %     [t,y]=ode15s(@funprke_exact,[time_beg time_end],X,options,t_ref,exact_arr);
% %     
% %     X=y(end,:);
% %     p(it+1) = X(1);
% % 
% % end


    time_beg = 0;
    time_end = dt*ntimes;
    X_beg = X;
    
    % tolerances for odesolvers
    rtol = 3e-14; abso = 3e-14;
    atol  = abso*ones(length(X),1);
    options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);    
    
    % call ode solver with prke function that linearly interpolates
    [t,y]=ode15s(@funprke_exact,[time_beg time_end],X,options,t_ref,exact_arr);
    
    X=y(end,:);
    p = X(1);
