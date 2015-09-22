function [X,dpdt,t,pow] =  solve_prke_ode(X,dt_macro,time_end,shape_beg,shape_end)

% make the problem-data a global variable
global dat

% create a list of times when reactivity will be recomputed for this macro
% time step
time_beg = time_end - dt_macro;

dat.ode.time_beg = time_beg;
dat.ode.time_end = time_end;
dat.ode.shape_beg = shape_beg;
dat.ode.shape_end = shape_end;

% tolerances for odesolvers
rtol = 3e-14; abso = 3e-14;
atol  = abso*ones(length(X),1);
options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);

% ODE solve for PRKEs
[t,y]=ode15s(@funprke,[time_beg time_end],X,options);
X=(y(end,:))';
pow=y(:,1);

% save value for dpdt at the end of the macro time step
dXdt=funprke(time_end,X);
dpdt=dXdt(1);

return
end


