function [ u ] = solve_TD_diffusion(u,dt,time_end)

global npar

TR = assemble_transient_operator(time_end);
M  = assemble_time_dependent_operator(time_end);

% build rhs from backward Euler time discretization
rhs = M*u;
% build system matrix
A = M-dt*TR;
if npar.set_bc_last
    [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
else
    rhs=apply_BC_vec_only(rhs);
end
% solve: M(unew-uold)/dt=TR.unew
u = A\rhs;


end

