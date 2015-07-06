function M = assemble_transient_operator_iqs(curr_time,p,dpdt)

global dat npar

n   = npar.n;
nnz = npar.nnz;
M   = sparse(2*n,2*n,2*(nnz+n));

D    = assemble_stiffness(dat.cdiff   ,curr_time);
A    = assemble_mass(     dat.siga    ,curr_time);
NFIp = assemble_mass(     dat.nusigf_p,curr_time) / npar.keff;
NFId = assemble_mass(     dat.nusigf_d,curr_time) / npar.keff * p;
Aiqs = assemble_mass(     dat.inv_vel ,curr_time) * dpdt/p;

% flux-flux matrix
tmp=NFIp-(D+A+Aiqs);
% prec-flux matrix
if ~npar.set_bc_last
    NFId=apply_BC_mat_only(NFId,npar.add_zero_on_diagonal);
    tmp =apply_BC_mat_only(tmp ,npar.add_zero_on_diagonal);
end

% generic matrix for precursors (it is diagonal because the unknown
% precursor concentrations are not the FEM expansion values but int b_i C
L=dat.lambda*speye(n);

% prec-prec matrix
Ldiag=L;
% flux-prec matrix
Loffd=L/p;

if ~npar.set_bc_last
    Ldiag=apply_BC_mat_only(L,npar.add_zero_on_diagonal);
    Loffd=apply_BC_mat_only(L,npar.add_zero_on_diagonal);
end

% finally, build the rhs operator for the transient problem
M(1:n    ,1:n    ) = tmp;
M(1:n    ,n+1:2*n) =  Loffd;
M(n+1:2*n,1:n    ) = NFId;
M(n+1:2*n,n+1:2*n) = -Ldiag;

return
end