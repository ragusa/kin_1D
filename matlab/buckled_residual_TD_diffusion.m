function F = buckled_residual_TD_diffusion(curr_time,u)

global npar
sav_value = npar.set_bc_last;
npar.set_bc_last = false;

phi = u(1:npar.n);
C = u(npar.n+1:npar.n*2);
T = u(npar.n*2+1:end);

TR = assemble_buckled_transient_operator(curr_time,T);
M  = assemble_buckled_time_dependent_operator(curr_time);

% build SS residual 
F = M\(TR*u);
% apply bc
F=apply_BC_vec_only(F);

npar.set_bc_last = sav_value;

end

