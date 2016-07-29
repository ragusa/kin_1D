function [X,dpdt,t,y,T] =  solve_prke_iqs_buckled(X,dt_macro,time_end,shape_beg,shape_end,Told,n_micro,freq_react)

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

% storage for t and y (for post-processing)
t(1) = time_beg;
y(1) = X(1);
etol = 1e-12;

% loop over micro time steps
for i=1:n_react
    t_beg = times_react_update(i);
    t_end = times_react_update(i+1);
    % weights for linear interpolation of the shape
    w1 = (time_end-t_beg)/dt_macro;
    w2 = (t_beg-time_beg)/dt_macro;
    % shape interpolation
    shape = shape_beg * w1 + shape_end * w2;
    X_old = X;
    t = [t t_end];
    y = [y X(1)];
    err = etol+1;
    mp_iter = 0;
    while err > etol
        mp_iter = mp_iter+1;
        X_prev = X;
        X = X_old;
        T = compute_temp_IQS(Told,[shape_beg shape],t,y);
        [rho_MGT(i+1),beff_MGT(i+1)]=compute_buckled_prke_parameters(t_end,shape,T);
        rho_MGT_iqs  = linspace(rho_MGT(i) ,rho_MGT(i+1) ,freq_react+1);
        beff_MGT_iqs = linspace(beff_MGT(i),beff_MGT(i+1),freq_react+1);
        for it=1:freq_react
            % build prke matrix
            A=[(rho_MGT_iqs(it+1)-beff_MGT_iqs(it+1)) dat.lambda ; ...
                beff_MGT_iqs(it+1)                   -dat.lambda];
            % solve
            X=(eye(2)-dt*A)\X;
            % storage
            t((i-1)*freq_react+it+1)=t((i-1)*freq_react+it)+dt;
            y((i-1)*freq_react+it+1)=X(1);
        end
        err = abs((X(1) - X_prev(1))/X_prev(1));
    end
%     fprintf('  number of temperature iterations = %d \n',mp_iter);
end
% fprintf('theta = %g, power= %g \n',theta,X(1));

% save value for dpdt at the end of the macro time step
dXdt = A*X;
dpdt=dXdt(1);

