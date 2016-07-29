function M = assemble_buckled_time_dependent_operator(curr_time)

global dat npar

n   = npar.n;
nnz = npar.nnz;
M   = sparse(3*n,3*n,nnz+2*n);

% flux-flux
IV = assemble_mass(dat.inv_vel,curr_time);
if ~npar.set_bc_last
    IV = apply_BC_mat_only(IV,npar.add_ones_on_diagonal);
end

% prec-prec
I=speye(2*n);

% save tr2.mat IV I;

% assemble time dep operator
M(1:n    ,1:n    ) = IV;
M(n+1:3*n,n+1:3*n) = I;

return
end