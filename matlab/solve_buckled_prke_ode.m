function [X,t,pow,T] =  solve_buckled_prke_ode(X,dt_macro,time_end,shape_beg,shape_end,Told,n_react)

% make the problem-data a global variable
global dat npar

% create a list of times when reactivity will be recomputed for this macro
% time step
time_beg = time_end - dt_macro;
times_react_update = time_beg + linspace(0,n_react,n_react+1) * dt_macro/n_react;

% tolerances for odesolvers
rtol = 3e-14; abso = 3e-14;
atol  = abso*ones(length(X),1);
options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);

etol = 1e-6;
[dat.ode.rho_MGT_beg,dat.ode.beff_MGT_beg]=compute_buckled_prke_parameters(time_beg,shape_beg,Told);
told = time_beg;
yold = X(1);

for i=1:npar.n_react
    dat.ode.time_beg = times_react_update(i);
    dat.ode.time_end = times_react_update(i+1);
    t = dat.ode.time_end;
    % weights for linear interpolation of the shape
    w1 = (time_end-dat.ode.time_end)/dt_macro;
    w2 = (dat.ode.time_end-time_beg)/dt_macro;
    % shape interpolation
    shape = shape_beg * w1 + shape_end * w2;
    X_old = X;
    y = X(1);
    err = etol+1;
    mp_iter = 0;
    while err > etol
        mp_iter = mp_iter+1;
        X_prev = X;
        X = X_old;
        T = compute_temp_IQS(Told,[shape_beg shape],[told; t],[yold; y(:,1)]);
        [dat.ode.rho_MGT_end,dat.ode.beff_MGT_end]=compute_buckled_prke_parameters(dat.ode.time_end,shape,T);
        % call ode solver with prke function that linearly interpolates
        [t,y]=ode15s(@funprke_lin_interp,[dat.ode.time_beg dat.ode.time_end],X,options);
        t = t(2:end);
        y = y(2:end,:);
        X=(y(end,:))';
        err = abs((X(1) - X_prev(1))/X_prev(1));
%     fprintf('  number of temperature iterations = %d \n',mp_iter);
    end
    told = [told; t]; t = [];
    yold = [yold; y(:,1)]; y = [];
    dat.ode.rho_MGT_beg = dat.ode.rho_MGT_end;
    dat.ode.beff_MGT_beg = dat.ode.beff_MGT_end;
end
% save data
t = told;
pow=yold;

% % save value for dpdt at the end of the macro time step
% if npar.solve_prke_compute_rho_each_time
%     dXdt=funprke(time_end,X);
% else
%     % react coef have already been computed. no need to re-compute them in
%     % the call to the prke function
%     switch npar.iqs_prke_interpolation_method
%         case 1
%             dXdt=funprke_lin_interp(time_end,X);
%         case 2
%             dXdt=funprke_linlin_interp(time_end,X);
%         case 3
%             dXdt=funprke_linhermite_interp(time_end,X);
%     end
% end
% dpdt=dXdt(1);

return
end


