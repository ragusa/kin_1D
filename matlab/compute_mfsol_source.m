function [S_phi_MGT,S_C_MGT]=compute_mfsol_source(curr_time,shape)

global dat npar
phi_adjoint = npar.phi_adj;

S_phi_mat = assemble_source( @mf_Source ,curr_time, true);
S_C_mat   = assemble_source( @mf_Source ,curr_time, false);
IV    = assemble_mass(dat.inv_vel ,curr_time);


IV   = apply_BC_mat_only(IV,npar.add_zero_on_diagonal);

S_phi = phi_adjoint' * S_phi_mat;
S_C = phi_adjoint' * S_C_mat;
MGT  = phi_adjoint' * IV * shape; % obviously, MGT = npar.K0 as well


S_phi_MGT  = S_phi/MGT;
S_C_MGT = S_C/MGT;


% prec_MGT = phi_adjoint' * C / MGT; 
% % obviously, this is the same as beff_MGT/dat.lambda because, at steady
% % state, C = NFid * phi0 / lambda, thus: 
% % phi_adjoint' * C = beff / lambda, and finally
% % phi_adjoint' * C / MGT = beff_MGT /lambda 


return
end