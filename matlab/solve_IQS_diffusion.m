function [u_shape, X, IV, err] = solve_IQS_diffusion(u_shape,X,dt_macro,time_end,n_micro,freq_react,theta,theta_log)

global dat npar

max_iter_iqs = 1;
tol_iqs=1e-7;

IV = assemble_mass(dat.inv_vel,time_end);

% save values at beginning of macro time step: they are needed in the IQS iteration 
X_beg=X;
u_shape_beg=u_shape;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
shape_beg=u_shape(1:npar.n);
shape_end=shape_beg;

for iter = 1: max_iter_iqs
    
    if npar.mf_sol
        S = assemble_Source_vector(time_end);
    end
    % solve for amplitude function
    [X,dpdt] =  solve_prke_iqs(X_beg,dt_macro,time_end,shape_beg,shape_end,n_micro,freq_react,theta,theta_log);
%     [X(1) dpdt]
    % assemble IQS    
    TR = assemble_transient_operator_iqs(time_end,X(1),dpdt);
    M  = assemble_time_dependent_operator(time_end);
     
    % build rhs from backward Euler time discretization
    rhs = M*u_shape_beg;
    % build system matrix
    A = M-dt_macro*TR;
    if npar.mf_sol
        rhs = rhs + S*dt_macro;
    end
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
    fprintf('  IQS iter %d, err %g, dlogp_dt %g \n',iter,err, dpdt/X(1));
    if err<tol_iqs
        break
    elseif ~theta_log
        u_shape = u_shape / ((npar.phi_adj)'*IV*shape_end/npar.K0);
        shape_end=u_shape(1:npar.n);
    end
end
% dat.dlogp_dt = [dat.dlogp_dt; dpdt/X(1)];

if err>=tol_iqs
    warning('IQS did not converge');
end

