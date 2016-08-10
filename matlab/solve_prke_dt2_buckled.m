function [X,dpdt,t,y,T] =  solve_prke_dt2_buckled(X,dt_macro,time_end,shape_beg,shape_end,Told,n_react)

% make the problem-data a global variable
global dat npar
rk = compute_SDIRKparams(3);

% number of reactivity updates during one macro time step
if abs(floor(n_react)-n_react) > 1e-10
    error('n_micro / freq_react must be an integer');
end
% create a list of times when reactivity will be recomputed for this macro
% time step
time_beg = time_end - dt_macro;
times_react_update = time_beg + linspace(0,n_react,n_react+1) * dt_macro/n_react;

[rho_MGT(1),beff_MGT(1)]=compute_buckled_prke_parameters(time_beg,shape_beg,Told);

% storage for t and y (for post-processing)
told = time_beg;
yold = X(1);
etol = 1e-8;
etol_dt2 = 1e-10;
dt_max = 10;

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
    t = t_end;
    y = X(1);
    err = etol+1;
    mp_iter = 0;
    while err > etol
        mp_iter = mp_iter+1;
        X_prev = X;
        X = X_old;
        T = compute_temp_IQS(Told,[shape_beg shape],[told t],[yold y]);
        [rho_MGT(i+1),beff_MGT(i+1)]=compute_buckled_prke_parameters(t_end,shape,T);
        rho_MGT_iqs   = @(t_cur) rho_MGT(i+1) *(t_cur-t_beg)/(t_end-t_beg) + rho_MGT(i) *(t_end-t_cur)/(t_end-t_beg);
        beff_MGT_iqs  = @(t_cur) beff_MGT(i+1)*(t_cur-t_beg)/(t_end-t_beg) + beff_MGT(i)*(t_end-t_cur)/(t_end-t_beg);
        A = @(tt) [(rho_MGT_iqs(tt)-beff_MGT_iqs(tt))  dat.lambda ; ...
                    beff_MGT_iqs(tt)                  -dat.lambda];
        t = t_beg;
        y = X(1);
        dt = (t_end-t_beg)*1e-4;
        keep_going=true;
        while (keep_going)
            t(end+1) = t(end) + dt;
            err_dt2 = etol_dt2+1;
            while err_dt2 > etol_dt2
                t(end) = t(end-1) + dt;
                dt_old = dt;
%                 % build prke matrix
%                 A_half=[(rho_MGT_iqs(mean(t(end-1:end)))-beff_MGT_iqs(mean(t(end-1:end))))  dat.lambda ; ...
%                                      beff_MGT_iqs(mean(t(end-1:end)))                      -dat.lambda];
%                 A=[(rho_MGT_iqs(t(end))-beff_MGT_iqs(t(end)))  dat.lambda ; ...
%                     beff_MGT_iqs(t(end))                      -dat.lambda];
%                 % solve
%                 X_half=(eye(2)-dt/2*A_half)\X;
%                 X_half=(eye(2)-dt/2*A)\X_half;
%                 X_full=(eye(2)-dt*A)\X;
                X_half = SDIRK33_solve(X,A,(t(end)+t(end-1))/2,dt/2,rk);
                X_half = SDIRK33_solve(X_half,A,t(end),dt/2,rk);
                X_full = SDIRK33_solve(X,A,t(end),dt,rk);
                err_dt2 = abs(X_half(1)-X_full(1))/X_half(1);
                dt = dt*(etol_dt2/err_dt2)^(1/4);
                if dt/dt_old > dt_max, dt=dt_old*dt_max; end
                dt;
            end
            X = X_half;
            if t(end)>=t_end, keep_going=false;
            elseif t(end)+dt>=t_end, dt = t_end - t(end); end
            % storage
            y(end+1)=X(1);
        end
        t = t(2:end);
        y = y(2:end);
        err = abs((X(1)-X_prev(1))/X_prev(1));
    end
    told = [told t]; t=[];
    yold = [yold y]; y=[];
%     fprintf('  number of temperature iterations = %d \n',mp_iter);
end
% fprintf('theta = %g, power= %g \n',theta,X(1));
t = told';
y = yold';

% save value for dpdt at the end of the macro time step
dXdt = A(t(end))*X;
dpdt=dXdt(1);

function X = SDIRK33_solve(X,A,time_end,dt,rk)
var = length(X);
k=zeros(var,rk.s);
Y=zeros(var,rk.s);
tn = time_end - dt;
ti = tn+rk.c*dt;

for n=1:rk.s
    M = eye(var)-rk.a(n,n)*dt*A(ti(n));
    rhs = X;
    for j=1:n-1,
        rhs = rhs + dt*rk.a(n,j)*k(:,j);
    end
    Y(:,n) = M\rhs;
    k(:,n) = A(ti(n))*Y(:,n);
end
X=Y(:,end);

end

end