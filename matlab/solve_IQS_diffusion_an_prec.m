function [u_shape, X,t,y] = solve_IQS_diffusion_an_prec(u_shape,X,dt_macro,time_end)

global dat npar

% shortcuts
lambda = dat.lambda;
dt = dt_macro;
C_old = u_shape(npar.n+1:end);

max_iter_iqs = 10;
tol_iqs=1e-7;
npar.theta_old=[];

% save values at beginning of macro time step: they are needed in the IQS iteration
X_beg=X;
u_shape_beg=u_shape;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
shape_beg=u_shape(1:npar.n);
shape_end=shape_beg;

for iter = 1: max_iter_iqs
    
    % solve for amplitude function
    if strcmpi(npar.prke_solve,'matlab')
        [X,dpdt,t,y] =  solve_prke_ode(X_beg,dt_macro,time_end,shape_beg,shape_end);
    else
        [X,dpdt,t,y] =  solve_prke_iqs(X_beg,dt_macro,time_end,shape_beg,shape_end,n_micro,freq_react);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assemble IQS
    
    % shortcut
    p=X(1);
    % interpolating polynomial
    if strcmpi(npar.prke_solve,'matlab')
        pp = interp1(t,y,'pchip','pp');
    else
        pp = interp1(t,y,'linear','pp');
    end
    f = @(t) ppval(pp,t);
    
    TR  = sparse(npar.n,npar.n,npar.nnz);
    
    D    = assemble_stiffness(dat.cdiff   ,time_end);
    A    = assemble_mass(     dat.siga    ,time_end);
    NFIp = assemble_mass(     dat.nusigf_p,time_end) / npar.keff;
    IV   = assemble_mass(     dat.inv_vel ,time_end);
    Aiqs = assemble_mass(     dat.inv_vel ,time_end) * dpdt/p;
    
    NFId_old = assemble_mass(dat.nusigf_d,time_end-dt) / npar.keff ;
    NFId_new = assemble_mass(dat.nusigf_d,time_end)    / npar.keff ;
    
    % flux-flux matrix
    TR=NFIp-(D+A+Aiqs);
    
    % integrals for the analytical expressions for the precursors
    t2=time_end;
    t1=t2-dt;
    
    A1= @(t)( ((t2-t)/dt).^2      .*exp(-lambda*(t2-t)) .*f(t) );
    A2= @(t)( (t-t1).*(t2-t)/dt^2 .*exp(-lambda*(t2-t)) .*f(t) );
    A3= @(t)( ((t-t1)/dt).^2      .*exp(-lambda*(t2-t)) .*f(t) );
    
    a1= quad(@(t)A1(t),t1,t2);
    a2= quad(@(t)A2(t),t1,t2);
    a3= quad(@(t)A3(t),t1,t2);
    
    
    % build transient matrix (divide lambda by p)
    TR = TR + lambda/p * ( a2*NFId_old + a3*NFId_new ) ;
    
    % build rhs (divide lambda by p)
    rhs = IV*shape_beg + dt * lambda/p *( C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*shape_beg );
    
    % build system matrix
    A = IV-dt*TR;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if npar.set_bc_last
        [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
    else
        rhs=apply_BC_vec_only(rhs);
    end
    % solve for new shape_end: M(unew-uold)/dt=TR.unew
    shape_end = A\rhs;
    
    % update precursors
    C_new =  C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*shape_beg + ( a2*NFId_old + a3*NFId_new )*shape_end ;
    % re-package as single solution vector
    u_shape = [ shape_end ; C_new];
    
    % check for tolerance
    err = abs( (npar.phi_adj)'*IV*shape_end/npar.K0  - 1);
    fprintf('  IQS iter %d, err %g \n',iter,err);
    if err<tol_iqs
        break
    else
        %         u_shape = u_shape / ((npar.phi_adj)'*IV*shape_end/npar.K0);
        %         shape_end=u_shape(1:npar.n);
    end
end

if err>=tol_iqs
    warning('IQS did not converge');
end

