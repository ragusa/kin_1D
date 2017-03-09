function F = residual_TD_diffusion(curr_time,u)

global npar
sav_value = npar.set_bc_last;
npar.set_bc_last = false;

TR = assemble_transient_operator(curr_time);
M  = assemble_time_dependent_operator(curr_time);

% build SS residual 
F = M\(TR*u);
% apply bc
F=apply_BC_vec_only(F);

npar.set_bc_last = sav_value;

end

