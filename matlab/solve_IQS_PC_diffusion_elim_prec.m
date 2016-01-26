function [u_end, X,t,y] = solve_IQS_PC_diffusion_elim_prec(u,X,dt_macro,time_end)

global io dat npar

% shortcuts
lambda = dat.lambda;
C_old = u(npar.n+1:end);

max_iter_iqs = npar.max_iter_iqs;
tol_iqs      = npar.tol_iqs;

npar.theta_old=[];

% save values at beginning of macro time step: they are needed in the IQS iteration
X_beg=X;
u_beg=u;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
% IV   = assemble_mass(     dat.inv_vel ,time_end);
z = (npar.phi_adj)'*npar.IV*u(1:npar.n)/npar.K0;
shape_beg=u(1:npar.n)/z;

for iter = 1: max_iter_iqs
    
    % solve time-dependent diffusion for flux
    [u_end, auxphi] = solve_TD_diffusion_elim_prec(u_beg,dt_macro,time_end);

    % get a predicted value of the end flux
    flux_end = u_end(1:npar.n);
    C_new_tmp= u_end(npar.n+1:end);
    
    % get a shape
    z = (npar.phi_adj)'*npar.IV*flux_end/npar.K0;
    shape_end = u_end(1:npar.n) / z;
    
    % solve for amplitude function
    if strcmpi(npar.prke_solve,'matlab')
        [X,w,t,y] =  solve_prke_ode(X_beg,dt_macro,time_end,shape_beg,shape_end);
%         [X,dpdt,t,y] =  solve_prke_ode_events(X_beg,dt_macro,time_end,shape_beg,shape_end);
    else
        [X,dpdt,t,y] =  solve_prke_iqs(X_beg,dt_macro,time_end,shape_beg,shape_end);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assemble IQS
    
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
%     % extract details from piece-wise polynomial by breaking it apart
%     [breaks,coefs,l,k,d] = unmkpp(pp);
%     % make the polynomial that describes the derivative
%     pp_der = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
    
    
    p = @(t) ppval(pp,t);   
    
%     [X_beg(1) X(1) y(end,1) p(time_end)]
    
    % integrals for the analytical expressions for the precursors
    t2=time_end;
    dt = dt_macro;
    t1=t2-dt;
    
    if npar.iqs_prke_interpolation_method<3
        A1= @(t)( ((t2-t)/dt).^2      .*exp(-lambda*(t2-t)) .*p(t) );
        A2= @(t)( (t-t1).*(t2-t)/dt^2 .*exp(-lambda*(t2-t)) .*p(t) );
        A3= @(t)( ((t-t1)/dt).^2      .*exp(-lambda*(t2-t)) .*p(t) );
        
        a1= integral(@(t)A1(t),t1,t2,'Reltol',eps);
        a2= integral(@(t)A2(t),t1,t2,'Reltol',eps);
        a3= integral(@(t)A3(t),t1,t2,'Reltol',eps);
        
        NFId_old = assemble_mass(dat.nusigf_d,time_end-dt) / npar.keff ;
        NFId_new = assemble_mass(dat.nusigf_d,time_end)    / npar.keff ;
        
        % update precursors
        C_new =  C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*shape_beg + ( a2*NFId_old + a3*NFId_new )*shape_end ;
        
    else
%                 for k=1:4
%                     Abeg{k}= @(t)( (t2-t)/dt.*t.^(4-k).*exp(-lambda*(t2-t)) .*p(t) );
%                     Aend{k}= @(t)( (t-t1)/dt.*t.^(4-k).*exp(-lambda*(t2-t)) .*p(t) );
%                     abeg(k)= integral(@(t)Abeg{k}(t),t1,t2,'Reltol',eps);
%                     aend(k)= integral(@(t)Aend{k}(t),t1,t2,'Reltol',eps);
%                 end
%                 NFId_old = assemble_mass(dat.nusigf_d,time_end-dt) / npar.keff ;
%                 NFId_new = assemble_mass(dat.nusigf_d,time_end)    / npar.keff ;
%         
%                 % update precursors
%                 C_new =  C_old*exp(-lambda*dt);
%                 for k=1:4
%                     C_new = C_new + (abeg(k)*NFId_old + aend(k)*NFId_new)*w(:,k);
%                 end
        time_beg=time_end-dt;
        rk=npar.rk;
        for i=1:rk.s
            ti = time_beg + rk.c(i)*dt;
            NFId_new = assemble_mass( dat.nusigf_d,ti) / npar.keff;
%             NFId_new=apply_BC_mat_only(NFId_new,npar.add_ones_on_diagonal);
            flux_ti = auxphi(1:npar.n,i);
            % get a shape
            z = (npar.phi_adj)'*npar.IV*flux_ti/npar.K0;
            shape_ti = flux_ti / z;
            % get better flux
            flux_ti = shape_ti * p(ti);
            % prec
            Ci = C_old;
            for j=1:i-1
                Ci = Ci + rk.a(i,j)*dt*fC(:,j);
            end
            Ci = Ci + rk.a(i,i)*dt * NFId_new*flux_ti;
            deno = (1+lambda*rk.a(i,i)*dt);
            Ci = Ci/deno;
            % store for temp SDIRK quantities
            if i<rk.s
                fC(:,i)=-lambda*Ci + NFId_new*flux_ti;
            end
        end
        C_new=Ci;
    end
    
    % re-package as single solution vector
    u_end = [ X(1)*shape_end ; C_new];
%     u_end = [ X(1)*shape_end ; C_new_tmp];
    u_end = [flux_end ; C_new_tmp];
%     [X(1) z abs(X(1)-z) abs(X(1)-p(time_end)) norm(C_new-C_new_tmp)]
    % check for tolerance
    err = abs( ((npar.phi_adj)'*npar.IV*shape_end)/npar.K0  - 1);
    if io.console_print
        fprintf('  IQS iter %d, err %g \n',iter,err);
    end
    if err<tol_iqs
        break
    else
        %         u_shape = u_shape / ((npar.phi_adj)'*IV*shape_end/npar.K0);
        %         shape_end=u_shape(1:npar.n);
    end
end

if err>=tol_iqs
    warning('IQS did not converge in %s',mfilename);
end