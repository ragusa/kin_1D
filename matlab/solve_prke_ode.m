function [X,dpdt,t,pow] =  solve_prke_ode(X,dt_macro,time_end,shape_beg,shape_end)

% make the problem-data a global variable
global dat npar

% create a list of times when reactivity will be recomputed for this macro
% time step
time_beg = time_end - dt_macro;

dat.ode.time_beg = time_beg;
dat.ode.time_end = time_end;

% tolerances for odesolvers
rtol = 3e-14; abso = 3e-14;
atol  = abso*ones(length(X),1);
options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);

% ODE solve for PRKEs
if npar.solve_prke_compute_rho_each_time
    % save shape to retrieve them in funprke
    dat.ode.shape_beg = shape_beg;
    dat.ode.shape_end = shape_end;
    % call ode solver with prke function that re-compute rho/beta at each time (expensive)
    [t,y]=ode15s(@funprke,[time_beg time_end],X,options);
else
    % compute prke parameters at beg/end time and linearly interpolate
    [dat.ode.rho_MGT_beg,dat.ode.beff_MGT_beg]=compute_prke_parameters(time_beg,shape_beg);
    [dat.ode.rho_MGT_end,dat.ode.beff_MGT_end]=compute_prke_parameters(time_end,shape_end);
    % call ode solver with prke function that linearly interpolates
    [t,y]=ode15s(@funprke_lin_interp,[time_beg time_end],X,options);
end
% save data
X=(y(end,:))';
pow=y(:,1);

% save value for dpdt at the end of the macro time step
if npar.solve_prke_compute_rho_each_time
    dXdt=funprke(time_end,X);
else
    % react coef have already been compute. no need to re-compute them in
    % the call to the prke function
    dXdt=funprke_lin_interp(time_end,X);
end
dpdt=dXdt(1);

return
end


