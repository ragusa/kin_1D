function [u_shape, X] = solve_IQS_diffusion(u_shape,X,dt_macro,time_end,n_micro,freq_react)

global dat npar

max_iter_iqs = 10;
tol_iqs=1e-5;

IV = assemble_mass(dat.inv_vel,time_end);

shape_beg=u_shape(1:npar.n);
shape_end=shape_beg;

for iter = 1: max_iter_iqs
    
    % solve for amplitude function
    [X,dpdt] =  solve_prke_iqs(X,dt_macro,time_end,shape_beg,shape_end,n_micro,freq_react);
    
    % assemble IQS    
    TR = assemble_transient_operator_iqs(time_end,X(1),dpdt);
    M  = assemble_time_dependent_operator(time_end);
    % build rhs from backward Euler time discretization
    rhs = M*u_shape;
    % build system matrix
    A = M-dt_macro*TR;
    if npar.set_bc_last
        [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
    else
        rhs=apply_BC_vec_only(rhs);
    end
    % solve for new shape_end: M(unew-uold)/dt=TR.unew
    u_shape = A\rhs;
    shape_end=u_shape(1:npar.n);
    
    % check for tolerance 
    err = abs( (npar.phi_adj)'*IV*shape_end/npar.K0  - 1);
    fprintf('  IQS iter %d, err %g \n',iter,err);
    if err<tol_iqs
        break
    end
end

if err>=tol_iqs
    error('IQS did not converge');
end

