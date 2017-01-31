function dydt=funprke_lag_interp(time,y)

% global
global dat npar

tn = dat.ode.tn;

% reactivity coef.  interpolation
rho_MGT = 0.0;
beff_MGT = 0.0;
for i=1:length(tn)
    rho = dat.ode.rho_MGT(i);
    beff = dat.ode.beff_MGT(i);
    for j=1:length(tn)
        if i~=j
            rho = rho.*(time-tn(j))/(tn(i)-tn(j));
            beff = beff.*(time-tn(j))/(tn(i)-tn(j));
        end
    end
    rho_MGT = rho_MGT + rho;
    beff_MGT = beff_MGT + beff;
end

% compute PRKE matrix at current time
J=[(rho_MGT-beff_MGT)   dat.lambda ; ...
    beff_MGT           -dat.lambda ];

S    = assemble_source(   dat.source_phi,time);
q_MGT    = (npar.phi_adj' * S)/npar.K0;

q = [q_MGT; 0];

% provide dydt to ode solver
dydt=J*y + q;

return
end
