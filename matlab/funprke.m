function dydt=funprke(time,y)

% global
global dat

time_beg = dat.ode.time_beg;
time_end = dat.ode.time_end;
dt_macro = time_end - time_beg;

w1 = (time_end-time)/dt_macro;
w2 = (time-time_beg)/dt_macro;
% shape interpolation
shape = dat.ode.shape_beg * w1 + dat.ode.shape_end * w2;

% compute prke parameters at current time
[rho_MGT,beff_MGT,q_MGT]=compute_prke_parameters(time,shape);


% compute PRKE matrix at current time
J=[(rho_MGT-beff_MGT)   dat.lambda ; ...
    beff_MGT           -dat.lambda ];

% provide dydt to ode solver
dydt=J*y + [q_MGT; 0];

return
end

