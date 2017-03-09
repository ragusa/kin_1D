function [C]=kinetics_init(phi,curr_time)

% make the problem-data a global variable
global dat npar 

NFId = assemble_mass(dat.nusigf_d,curr_time) /npar.keff;
if ~npar.set_bc_last
    NFId = apply_BC_mat_only(NFId,npar.add_zero_on_diagonal);
end
C = NFId*phi/dat.lambda;

% verif
% M=assemble_transient_operator(curr_time);
% plot(M*[phi;C])

% adjoint flux
npar.phi_adj=phi;

% normalization constant
IV = assemble_mass(dat.inv_vel,curr_time);
% npar.IV0=IV;
K0 = npar.phi_adj' * IV * phi;
npar.K0 = K0;

% total power
dat.Ptot = compute_power(dat.nusigf,curr_time,phi);

if ~npar.set_bc_last
    % add precursors BC values
    if dat.bc.left.type==2
        dat.bc.left.C(2) = C(1);
    end
    if dat.bc.rite.type==2
        dat.bc.rite.C(2) = C(end);
    end
end

return
end


