function [X,dpdt] =  solve_prke_iqs(X,dt_macro,time_end,shape_beg,shape_end,n_micro,freq_react)

% make the problem-data a global variable
global dat 

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

% make piecewise linear value for rho/beff at each micro time steps;
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
    A=[(rho_MGT(it)-beff_MGT(it)) dat.lambda ; ...
        beff_MGT                 -dat.lambda];
    % solve
    X=(eye(2)-dt*A)\X;
end

% save value for dpdt at the end of the macro time step
dXdt = A*X;
dpdt=dXdt(1);

