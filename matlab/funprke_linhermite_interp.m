function dydt=funprke_linhermite_interp(time,y)

% global
global dat npar io

nw = dat.ode.nw;

time_beg = dat.ode.time_beg;
time_end = dat.ode.time_end;
dt_macro = time_end - time_beg;

w1 = (time_end-time)/dt_macro;
w2 = (time-time_beg)/dt_macro;

% reactivity coef.  interpolation
rho_MGT  = 0;
beff_MGT = 0;
flux_time_t = 0;
for k=1:nw
    rho_MGT  = rho_MGT  + ( dat.ode.rho_MGT_beg(k)  * w1 + dat.ode.rho_MGT_end(k)  * w2 ) * time^(nw-k);
    beff_MGT = beff_MGT + ( dat.ode.beff_MGT_beg(k) * w1 + dat.ode.beff_MGT_beg(k) * w2 ) * time^(nw-k);
    flux_time_t = flux_time_t + io.w(:,k)*time^(nw-k);
end
% renromalize in case flux and not shape was used to compute prke parameters
renorm_flx = ( ((npar.phi_adj)'*npar.IV*flux_time_t) / npar.K0 )
rho_MGT  = rho_MGT  / renorm_flx ;
beff_MGT = beff_MGT / renorm_flx ;

% compute PRKE matrix at current time
J=[(rho_MGT-beff_MGT)   dat.lambda ; ...
    beff_MGT           -dat.lambda ];

% provide dydt to ode solver
dydt=J*y;

return
end

