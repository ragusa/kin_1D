function [rho_MGT,beff_MGT,prec_MGT]=compute_prke_parameters(curr_time,shape,C)

global dat npar

phi_adjoint = npar.phi_adj;

D    = assemble_stiffness(dat.cdiff   ,curr_time);
A    = assemble_mass(     dat.siga    ,curr_time);
NFI  = assemble_mass(     dat.nusigf  ,curr_time) / npar.keff;
NFId = assemble_mass(     dat.nusigf_d,curr_time) / npar.keff;
IV   = assemble_mass(     dat.inv_vel ,curr_time);

IV   = apply_BC_mat_only(IV,npar.add_zero_on_diagonal);
NFId = apply_BC_mat_only(NFId,npar.add_zero_on_diagonal);

PmM = apply_BC_mat_only(NFI-D-A,npar.add_zero_on_diagonal);
% rho  = phi_adjoint' * (NFI-D-A) * shape;
rho  = phi_adjoint' * (PmM) * shape;
beff = phi_adjoint' * NFId * shape;
MGT  = phi_adjoint' * IV * shape;

rho_MGT  = rho/MGT;
beff_MGT = beff/MGT;

prec_MGT = phi_adjoint' * C / MGT;


return
end