function [X,dpdt] =  solve_prke_iqs(X,dt_macro,time_end,shape_beg,shape_end,n_micro,freq_react,theta,theta_log)

% make the problem-data a global variable
global dat npar

% time step for prke solve
dt = dt_macro/n_micro;
% storage for power level at micro steps
pow_level=zeros(n_micro,1);

% number of reactivity updates during one macro time step
n_react = n_micro / freq_react;
if abs(floor(n_react)-n_react) > 1e-10
    error('n_micro / freq_react must be an integer');
end
% create a list of times when reactivity will be recomputed for this macro
% time step
time_beg = time_end - dt_macro;
times_react_update = time_beg + linspace(0,n_react,n_react+1) * dt_macro/n_react;
    if theta_log
        theta = ((npar.phi_adj)'*npar.IVel*(shape_end-shape_beg))/(dt*npar.K0);
        fprintf('theta = %g \n',theta);
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
    if npar.mf_sol
        [S_phi_MGT(i), S_C_MGT(i)]=compute_mfsol_source(times_react_update(i),shape);
    else
        S_phi_MGT(i) = 0;
        S_C_MGT(i)=0; 
    end
end

% make piecewise linear value for rho/beff at each micro time step;
% since we are using BE, we do not need the value at time_beg
rho_MGT_iqs  = [];
beff_MGT_iqs = [];
S_phi_MGT_iqs= [];
S_C_MGT_iqs  = [];
for i=1:n_react
    aux = linspace(rho_MGT(i),rho_MGT(i+1),freq_react+1);
    rho_MGT_iqs = [ rho_MGT_iqs aux(2:end)];
    aux = linspace(beff_MGT(i),beff_MGT(i+1),freq_react+1);
    beff_MGT_iqs = [ beff_MGT_iqs aux(2:end)];
    aux = linspace(S_phi_MGT(i),S_phi_MGT(i+1),freq_react+1);
    S_phi_MGT_iqs = [ S_phi_MGT_iqs aux(2:end)];
    aux = linspace(S_C_MGT(i),S_C_MGT(i+1),freq_react+1);
    S_C_MGT_iqs = [ S_C_MGT_iqs aux(2:end)];
end

% loop over micro time steps
for it=1:n_micro
    % build prke matrix
    A=[(rho_MGT_iqs(it)-beff_MGT_iqs(it)-theta), dat.lambda ; ...
        beff_MGT_iqs(it),                    -dat.lambda-theta];
    % rhs
    rhs = X + [S_phi_MGT_iqs(it) ; S_C_MGT_iqs(it)]*dt;
    % solve
    X=(eye(2)-dt*A)\rhs;
end

% save value for dpdt at the end of the macro time step
dXdt = A*X;
dpdt=dXdt(1);

