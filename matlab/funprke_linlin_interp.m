function dydt=funprke_linlin_interp(time,y)

% global
global dat

time_beg = dat.ode.time_beg;
time_end = dat.ode.time_end;
dt_macro = time_end - time_beg;

w1 = (time_end-time)/dt_macro;
w2 = (time-time_beg)/dt_macro;

% reactivity coef.  interpolation
rho_MGT  = ( dat.ode.rho_MGT_beg(1)  * w1 + dat.ode.rho_MGT_end(1)  * w2 ) * w1 + ...
          +( dat.ode.rho_MGT_beg(2)  * w1 + dat.ode.rho_MGT_end(2)  * w2 ) * w2;
beff_MGT = ( dat.ode.beff_MGT_beg(1) * w1 + dat.ode.beff_MGT_beg(1) * w2 ) * w1 + ...
          +( dat.ode.beff_MGT_beg(2) * w1 + dat.ode.beff_MGT_beg(2) * w2 ) * w2;

% compute PRKE matrix at current time
J=[(rho_MGT-beff_MGT)   dat.lambda ; ...
    beff_MGT           -dat.lambda ];

% provide dydt to ode solver
dydt=J*y;

return
end

