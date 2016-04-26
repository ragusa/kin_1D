function [u_shape, X,t,y] = solve_IQS_diffusion_elim_prec(u_shape,X,dt_macro,time_end)

global io dat npar

% shortcuts
lambda = dat.lambda;
dt = dt_macro;
C_old = u_shape(npar.n+1:end);

max_iter_iqs = npar.max_iter_iqs;
tol_iqs      = npar.tol_iqs;

npar.theta_old=[];
% save values at beginning of macro time step: they are needed in the IQS iteration
X_beg=X;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
shape_beg=u_shape(1:npar.n);
shape_end=shape_beg;

% SDIRK solve
% Yi = yn + h sum_j aij f(tj,Yj)
% gives (f=Qy+z)
%  [I-aii.h.Q(ti)]Yi = yn + h.aii zi + h sum(j=1..i-1) ( aij {Q(tj)Yj + zj} )
rk=npar.rk;
% beginning of the time interval
time_beg = time_end - dt;
% storage for temp SDIRK quantities (-lambda Ci + Bi Phii)
f=zeros(length(shape_beg),rk.s);
fC=f;
NFId_old = assemble_mass(dat.nusigf_d,time_beg) / npar.keff;


for iter = 1: max_iter_iqs
    
    % solve for amplitude function over the entire macro-step
    if strcmpi(npar.prke_solve,'matlab')
        [X,w,t,y] =  solve_prke_ode(X_beg,dt_macro,time_end,shape_beg,shape_end);
    else
        [X,dpdt,t,y] =  solve_prke_iqs(X_beg,dt_macro,time_end,shape_beg,shape_end,npar.n_micro,npar.freq_react);
    end
    
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
    
    % extract details from piece-wise polynomial by breaking it apart
    [breaks,coefs,l,k,d] = unmkpp(pp);
    % make the polynomial that describes the derivative
    pp_der = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:rk.s
        % compute stage time
        ti = time_beg + rk.c(i)*dt;
        
        % OLD shortcut
        % OLD p=X(1);
        pi = ppval(pp,ti);
        dpdti = ppval(pp_der,ti);
        
        % assemble IQS
        D    = assemble_stiffness(dat.cdiff   ,ti);
        A    = assemble_mass(     dat.siga    ,ti);
        NFIp = assemble_mass(     dat.nusigf_p,ti) / npar.keff;
        IV   = assemble_mass(     dat.inv_vel ,ti);
        Aiqs = assemble_mass(     dat.inv_vel ,ti) * dpdti/pi;
        NFId_new = assemble_mass( dat.nusigf_d,ti) / npar.keff;
        
        % flux-flux matrix
        TR=NFIp-(D+A+Aiqs);
        
        % expressions for the precursors
        % Ci = [ C_old + h aii FISdi.Phii + sum_{j<i} h aij(-lambda.Cj + FISdj.Phij)]/(1+lambda h aii)
        % let fci = -lambda.Ci + FISdi.Phii
        % in the above, Phii is the flux, not the shape since we are looking at the precursors eqs
        Ci = C_old;
        for j=1:i-1
            Ci = Ci + rk.a(i,j)*dt*fC(:,j);
        end
        deno = (1+lambda*rk.a(i,i)*dt);
        Ci = Ci/deno; % this is Ci without: h.aii.FISdi.Phii/deno
        
        % build transient matrix (divide lambda by p)
        TR = TR + (lambda * ( rk.a(i,i)*dt / deno )) * NFId_new;
        % build system matrix
        A = IV-rk.a(i,i)*dt*TR;
        
        % build rhs (divide lambda by p)
        zi = lambda/pi*Ci;
        rhs = IV*shape_beg + rk.a(i,i)*dt *zi;
        for j=1:i-1
            rhs = rhs + rk.a(i,j)*dt*f(:,j);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if npar.set_bc_last
            [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
        else
            rhs=apply_BC_vec_only(rhs);
        end
        % solve for new shape_end
        shape_end = A\rhs;
        % finish Ci
        Ci = Ci +  rk.a(i,i)*dt / deno * NFId_new*shape_end*pi;
        % store for temp SDIRK quantities
        f(:,i)=TR*shape_end + zi;
        fC(:,i)=-lambda*Ci + NFId_new*shape_end*pi;
    end % time integration loop end
    
    % save for hermite interp
    if npar.iqs_prke_interpolation_method>=3
        dat.ode.f_end=IV\f(:,rk.s);
    end
    
    % re-package as single solution vector
    u_shape = [ shape_end ; Ci];
    
    % check for tolerance
    err = abs( (npar.phi_adj)'*IV*shape_end/npar.K0  - 1);
    if io.console_print
        fprintf('  IQS iter %d, err %g \n',iter,err);
    end
    if err<tol_iqs
        break
    else
%                 u_shape = u_shape / ((npar.phi_adj)'*IV*shape_end/npar.K0);
%                 shape_end=u_shape(1:npar.n);
    end
end

if err>=tol_iqs
%     warning('IQS did not converge in %s',mfilename);
end


% renormalize anyway
u_shape = u_shape / ( ((npar.phi_adj)'*npar.IV*shape_end) / npar.K0 );