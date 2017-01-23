function [u_shape, X, t, y] = solve_IQS_diffusion(u_shape,X,dt_macro,tn)

global io npar dat

max_iter_iqs = npar.max_iter_iqs;
tol_iqs      = npar.tol_iqs;

npar.theta_old=[];

% IV = assemble_mass(dat.inv_vel,time_end);

% save values at beginning of macro time step: they are needed in the IQS iteration
X_beg=X(:,end);
u_shape_beg=u_shape;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
shape_beg=u_shape(1:npar.n,end);
shape_end=shape_beg;

% beginning of the time interval
time_end = tn(end);
time_beg = time_end - dt_macro;
rk=npar.rk;
if strcmp(npar.method,'BDF')
    order = length(tn);
    bdf = npar.bdf;
end
% storage for temp SDIRK quantities
K=zeros(length(u_shape(:,end)),rk.s);

for iter = 1: max_iter_iqs
    
    % solve for amplitude function
    if strcmpi(npar.prke_solve,'matlab')
        [X,dpdt,t,y] =  solve_prke_ode(X_beg,dt_macro,time_end,shape_beg,shape_end);
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
        ti = time_beg + rk.c(i)*dt_macro;
        
        pi = ppval(pp,ti);
        dpdti = ppval(pp_der,ti);
        
        % assemble IQS
        TR = assemble_transient_operator_iqs(ti,pi,dpdti);
        M  = assemble_time_dependent_operator(ti);
        S = assemble_source_operator(ti)/pi;
        
        if strcmp(npar.method,'BDF')
            % build system matrix
            A = M - dt_macro*bdf.b(order)*TR;
            % build rhs
            rhs = 0;
            for j=1:order
                rhs = rhs + bdf.a(order,j)*u_shape_beg(:,j);
            end
            rhs = M*rhs + dt_macro*bdf.b(order)*S;
        else
            % build system matrix
            A = M-rk.a(i,i)*dt_macro*TR;
            % build rhs 
            rhs = M*u_shape_beg + rk.a(i,i)*dt_macro*S;
            for j=1:rk.s-1
                rhs = rhs + rk.a(i,j)*dt_macro*K(:,j);
            end
        end
        % apply BC
        if npar.set_bc_last
            [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
        else
            rhs=apply_BC_vec_only(rhs);
        end
        % solve 
        u_shape = A\rhs;
        % store for temp SDIRK quantities
        TR=apply_BC_mat_only(TR,npar.add_ones_on_diagonal);
        K(:,i)=TR*u_shape + S;
    end % end rk loop
    % new shape
    shape_end=u_shape(1:npar.n);
    % save for hermite interp
    if npar.iqs_prke_interpolation_method>=3
        IV = assemble_mass(dat.inv_vel,time_end);
        dat.ode.f_end=IV\K(1:npar.n,rk.s);
    end
    
    % check for tolerance
    err = abs( ((npar.phi_adj)'*npar.IV*shape_end)/npar.K0  - 1);
    if io.console_print
        fprintf('  IQS iter %d, err %g \n',iter,err);
    end
    if err<tol_iqs
        break
    else
        %         u_shape = u_shape / ((npar.phi_adj)'*npar.IV*shape_end/npar.K0);
        %         shape_end=u_shape(1:npar.n);
    end
end

if err>=tol_iqs
%     warning('IQS did not converge in %s',mfilename);
end

% renormalize anyway
u_shape = u_shape / ( ((npar.phi_adj)'*npar.IV*shape_end) / npar.K0 );
