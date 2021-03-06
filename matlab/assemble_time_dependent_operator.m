function M = assemble_time_dependent_operator(curr_time)

global dat npar

n   = npar.n;
nnz = npar.nnz;
M   = sparse(2*n,2*n,nnz+n);

% flux-flux
IV = assemble_mass(dat.inv_vel,curr_time);
if ~npar.set_bc_last
    IV = apply_BC_mat_only(IV,npar.add_ones_on_diagonal);
end

% prec-prec
I=speye(n);

% save tr2.mat IV I;

% assemble time dep operator
M(1:n    ,1:n    ) = IV;
M(n+1:2*n,n+1:2*n) = I;

return
end