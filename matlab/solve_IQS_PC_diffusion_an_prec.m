function [u_end, X,t,y] = solve_IQS_PC_diffusion_an_prec(u,X,dt_macro,time_end)

global dat npar

% shortcuts
lambda = dat.lambda;
dt = dt_macro;
C_old = u(npar.n+1:end);

max_iter_iqs = npar.max_iter_iqs;
tol_iqs      = npar.tol_iqs;

npar.theta_old=[];

% save values at beginning of macro time step: they are needed in the IQS iteration
X_beg=X;
u_beg=u;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
IV   = assemble_mass(     dat.inv_vel ,time_end);
z = (npar.phi_adj)'*IV*u(1:npar.n)/npar.K0;
shape_beg=u(1:npar.n)/z;

for iter = 1: max_iter_iqs
    
    % solve time-dependent diffusion for flux
    u_end = solve_TD_diffusion_an_prec(u_beg,dt_macro,time_end);
    
    % get a predicted value of the end flux
    flux_end = u_end(1:npar.n);
    
    % get a shape
    z = (npar.phi_adj)'*IV*flux_end/npar.K0;
    shape_end = u_end(1:npar.n) / z;
    
    % solve for amplitude function
    if strcmpi(npar.prke_solve,'matlab')
        [X,dpdt,t,y] =  solve_prke_ode(X_beg,dt_macro,time_end,shape_beg,shape_end);
    else
        [X,dpdt,t,y] =  solve_prke_iqs(X_beg,dt_macro,time_end,shape_beg,shape_end);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assemble IQS
    
    % shortcut
    p=X(1);
    % interpolating polynomial
    if strcmpi(npar.prke_solve,'matlab')
        if npar.int_order==3
            pp = interp1(t,y,'spline','pp');
        elseif npar.int_order==2
            pp = interp1(t,y,'pchip','pp');
        elseif npar.int_order==1
            pp = interp1(t,y,'linear','pp');
        end        
    else
        pp = interp1(t,y,'linear','pp');
    end
    f = @(t) ppval(pp,t);   
    
    % integrals for the analytical expressions for the precursors
    t2=time_end;
    t1=t2-dt;
    
    A1= @(t)( ((t2-t)/dt).^2      .*exp(-lambda*(t2-t)) .*f(t) );
    A2= @(t)( (t-t1).*(t2-t)/dt^2 .*exp(-lambda*(t2-t)) .*f(t) );
    A3= @(t)( ((t-t1)/dt).^2      .*exp(-lambda*(t2-t)) .*f(t) );
    
    a1= integral(@(t)A1(t),t1,t2,'Reltol',eps);
    a2= integral(@(t)A2(t),t1,t2,'Reltol',eps);
    a3= integral(@(t)A3(t),t1,t2,'Reltol',eps);
    
    NFId_old = assemble_mass(dat.nusigf_d,time_end-dt) / npar.keff ;
    NFId_new = assemble_mass(dat.nusigf_d,time_end)    / npar.keff ;

    % update precursors
    C_new =  C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*shape_beg + ( a2*NFId_old + a3*NFId_new )*shape_end ;
    % re-package as single solution vector
    u_end = [ X(1)*shape_end ; C_new];
    
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
    warning('IQS_PC did not converge');
end

