function [u_shape, X,t,y] = solve_IQS_diffusion_elim_prec(u_shape,X,dt_macro,tn)

global io dat npar

% shortcuts
lambda = dat.lambda;
dt = dt_macro;
C_old = u_shape(npar.n+1:end,:);

max_iter_iqs = npar.max_iter_iqs;
tol_iqs      = npar.tol_iqs;

npar.theta_old=[];
% save values at beginning of macro time step: they are needed in the IQS iteration
X_beg=X;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
shape_beg=u_shape(1:npar.n,:);
shape_end=shape_beg(:,end);

% SDIRK solve
% Yi = yn + h sum_j aij f(tj,Yj)
% gives (f=Qy+z)
%  [I-aii.h.Q(ti)]Yi = yn + h.aii zi + h sum(j=1..i-1) ( aij {Q(tj)Yj + zj} )
rk=npar.rk;
if strcmp(npar.method,'BDF')
    bdf = npar.bdf;
    order = length(tn)-1;
end

% beginning of the time interval
time_end = tn(end);
time_beg = time_end - dt_macro;
% storage for temp SDIRK quantities (-lambda Ci + Bi Phii)
f=zeros(length(shape_beg),rk.s);
fC=f;
NFId_old = assemble_mass(dat.nusigf_d,time_beg) / npar.keff;


for iter = 1: max_iter_iqs
    shape_prev = shape_end;
    X_prev = X;
    % solve for amplitude function over the entire macro-step
    if strcmpi(npar.prke_solve,'matlab')
        [X,dXdt,w,t,y] =  solve_prke_ode(X_beg(:,end),dt_macro,tn,shape_beg(:,1:end),shape_end);
    else
        [X,dpdt,t,y] =  solve_prke_iqs(X_beg(:,end),dt_macro,time_end,shape_beg(:,end),shape_end,npar.n_micro,npar.freq_react);
    end
%     X(1) = (1+time_end);
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
%         pi = (1+ti);
%         dpdti = 1;
        pi = ppval(pp,ti);
%         dpdti = ppval(pp_der,ti);
%         tmp = funprke(ti,X);
%         dpdti = tmp(1);
        dpdti = dXdt(1);
        
        % assemble IQS
        D    = assemble_stiffness(dat.cdiff   ,ti);
        A    = assemble_mass(     dat.siga    ,ti);
        NFIp = assemble_mass(     dat.nusigf_p,ti) / npar.keff;
        IV   = assemble_mass(     dat.inv_vel ,ti);
        Aiqs = assemble_mass(     dat.inv_vel ,ti) * dpdti/pi;
        NFId_new = assemble_mass( dat.nusigf_d,ti) / npar.keff;
        S    = assemble_source(   dat.source_phi,ti)/pi;
        
        % flux-flux matrix
        TR=NFIp-(D+A+Aiqs);
        
        if strcmp(npar.method,'BDF')
            % expressions for the precursors
            Ci = 0;
            for j=1:order
                Ci = Ci + bdf.a(order,j)*C_old(:,j);
            end
            deno = (1 + dt*bdf.b(order)*lambda);
            Ci = Ci/deno;
            % build transient matrix
            TR = TR + dt*bdf.b(order)*lambda/deno*NFId_new;
            % build sysetm matrix
            A = IV - dt*bdf.b(order)*TR;
            % build rhs
            zi = lambda/pi*Ci+S;
            rhs = 0;
            for j=1:order
                rhs = rhs + bdf.a(order,j)*shape_beg(:,j);
            end
            rhs = dt*bdf.b(order)*zi + IV*rhs;
            
        else
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
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if npar.set_bc_last
            [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
        else
            rhs=apply_BC_vec_only(rhs);
        end
        % solve for new shape_end
        shape_end = A\rhs;
%         shape_end = npar.shape_exact(npar.x_dofs',ti);
%         shape_end = shape_beg(:,end);
        % finish Ci
        if strcmp(npar.method,'BDF')
            Ci = Ci +  bdf.b(order)*dt / deno * NFId_new*shape_end*pi;
        else
            Ci = Ci +  rk.a(i,i)   *dt / deno * NFId_new*shape_end*pi;
        end
        % store for temp SDIRK quantities
        f(:,i)=TR*shape_end + zi;
        fC(:,i)=-lambda*Ci + NFId_new*shape_end*pi;
    end % time integration loop end
    
    % save for hermite interp
    if npar.iqs_prke_interpolation_method>=3
%         dat.ode.f_end=IV\f(:,rk.s);
%         dat.ode.f_end=npar.dshape_exact(npar.x_dofs',time_end);
    end
    
    % re-package as single solution vector
    u_shape = [ shape_end ; Ci];
    
    % check for tolerance
    err = 0;
    if strcmp(npar.conv_criteria,'Linf')
        err = compute_conv_error(shape_end,shape_prev,Inf);
    elseif strcmp(npar.conv_criteria,'L2')
        err = compute_conv_error(shape_end,shape_prev,2);
    elseif max(strcmp(npar.conv_criteria,{'amp', 'rho', 'param_all'}))
        if max(strcmp(npar.conv_criteria,{'rho', 'param_all'}))
            [rho_end,~,~]=compute_prke_parameters(time_end,shape_end);
            [rho_prev,~,~]=compute_prke_parameters(time_end,shape_prev);
            err = max(err,abs(rho_end - rho_prev));
        end
        if max(strcmp(npar.conv_criteria,{'amp', 'param_all'}))
            err = max(err,abs(X(1) - X_prev(1)));
        end
        if strcmp(npar.conv_criteria,'param_all')
            err = max(err,abs( (npar.phi_adj)'*IV*shape_end/npar.K0  - 1));
        end
    else
        err = abs( (npar.phi_adj)'*IV*shape_end/npar.K0  - 1);
    end
    
    if io.console_print
        fprintf('  IQS iter %d, err %g \n',iter,err);
    end
    if err<tol_iqs
        break
    elseif npar.scale_IQS_each
        u_shape = u_shape / ((npar.phi_adj)'*IV*shape_end/npar.K0);
        shape_end=u_shape(1:npar.n);
    end
end

dat.iter_IQS = [dat.iter_IQS iter];
dat.iter_err = [dat.iter_err err];

if err>=tol_iqs
%     warning('IQS did not converge in %s',mfilename);
end


% renormalize anyway
if npar.scale_IQS
    u_shape = u_shape / ( ((npar.phi_adj)'*npar.IV*shape_end) / npar.K0 );
end