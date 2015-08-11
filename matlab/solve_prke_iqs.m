function [X,dpdt] =  solve_prke_iqs(X,dt_macro,time_end,shape_beg,shape_end,n_micro,freq_react)

% make the problem-data a global variable
global dat npar

% time step for prke solve
dt = dt_macro/n_micro;

% number of reactivity updates during one macro time step
n_react = n_micro / freq_react;
if abs(floor(n_react)-n_react) > 1e-10
    error('n_micro / freq_react must be an integer');
end
% create a list of times when reactivity will be recomputed for this macro
% time step
time_beg = time_end - dt_macro;
times_react_update = time_beg + linspace(0,n_react,n_react+1) * dt_macro/n_react;

theta_beg = npar.phi_adj' * npar.IV(1:npar.n,1:npar.n) * shape_beg;
theta_end = npar.phi_adj' * npar.IV(1:npar.n,1:npar.n) * shape_end;
theta = -(theta_end-theta_beg)/dt_macro/npar.K0;

if ~isempty(npar.theta_old)
    % damping
    a=1;
    theta = a*theta+(1-a)*npar.theta_old;
    theta=theta;
else
    npar.theta_old=theta;
end


% compute prke parameters at the above times 
for i=1:n_react+1 % same as length(times_react_update)
    % weights for linear interpolation of the shape
    w1 = (time_end-times_react_update(i))/dt_macro;
    w2 = (times_react_update(i)-time_beg)/dt_macro;
    % shape interpolation
    shape = shape_beg * w1 + shape_end * w2;
    % get prke values at different times
    [rho_MGT(i),beff_MGT(i)]=compute_prke_parameters(times_react_update(i),shape);
end

% make piecewise linear value for rho/beff at each micro time step;
% since we are using BE, we do not need the value at time_beg
rho_MGT_iqs  = [];
beff_MGT_iqs = [];
for i=1:n_react
    aux = linspace(rho_MGT(i),rho_MGT(i+1),freq_react+1);
    rho_MGT_iqs = [ rho_MGT_iqs aux(2:end)];
    aux = linspace(beff_MGT(i),beff_MGT(i+1),freq_react+1);
    beff_MGT_iqs = [ beff_MGT_iqs aux(2:end)];
end

% loop over micro time steps
for it=1:n_micro
    % build prke matrix
    A=[(rho_MGT_iqs(it)-beff_MGT_iqs(it)-theta) dat.lambda ; ...
        beff_MGT_iqs(it)                 -(dat.lambda+theta)];
    % solve
    X=(eye(2)-dt*A)\X;
end

fprintf('theta = %g, power= %g \n',theta,X(1));

% save value for dpdt at the end of the macro time step
dXdt = A*X;
dpdt=dXdt(1);

