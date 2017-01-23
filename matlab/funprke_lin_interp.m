function dydt=funprke_lin_interp(time,y)

% global
global dat

time_beg = dat.ode.time_beg;
time_end = dat.ode.time_end;
dt_macro = time_end - time_beg;

w1 = (time_end-time)/dt_macro;
w2 = (time-time_beg)/dt_macro;

% reactivity coef.  interpolation
rho_MGT  = dat.ode.rho_MGT_beg  * w1 + dat.ode.rho_MGT_end  * w2;
beff_MGT = dat.ode.beff_MGT_beg * w1 + dat.ode.beff_MGT_end * w2;
q_MGT    = dat.ode.q_MGT_beg    * w1 + dat.ode.q_MGT_end    * w2;


% compute PRKE matrix at current time
J=[(rho_MGT-beff_MGT)   dat.lambda ; ...
    beff_MGT           -dat.lambda ];

q = [q_MGT; 0];

% provide dydt to ode solver
dydt=J*y + q;

return
end

