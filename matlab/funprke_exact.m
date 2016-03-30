function dydt=funprke_exact(time,y,t_ref,exact_arr)

global dat npar

ind = find (t_ref<time);
if isempty(ind)
    i1=1;
else
    i1=ind(end);
end
ind = find (t_ref>=time);
if isempty(ind)
    i2=length(t_ref);
else
    i2=ind(1);
end
if( i1==i2)
    flux_time_t = exact_arr(1:npar.n,i1);
else
    dt = t_ref(i2) - t_ref(i1);
    w1 = (t_ref(i2)-time)/dt;
    w2 = (time-t_ref(i1))/dt;
    flux_time_t = w1*exact_arr(1:npar.n,i1) + w2*exact_arr(1:npar.n,i2);
end
renorm_flx = ( ((npar.phi_adj)'*npar.IV*flux_time_t) / npar.K0 );
shape = flux_time_t / renorm_flx;
[rho_MGT,beff_MGT] = compute_prke_parameters(time,shape);
 


% compute PRKE matrix at current time
J=[(rho_MGT-beff_MGT)   dat.lambda ; ...
    beff_MGT           -dat.lambda ];

% provide dydt to ode solver
dydt=J*y;

return
end

