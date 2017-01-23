function [u_shape, X,t,y] = solve_IQS_diffusion_an_prec(u_shape,X,dt_macro,tn)

global io dat npar

% shortcuts
lambda = dat.lambda;
dt = dt_macro;
C_old = u_shape(npar.n+1:end,end);
rk=npar.rk;
if strcmp(npar.method,'BDF')
    bdf = npar.bdf;
end
order = length(tn);

max_iter_iqs = npar.max_iter_iqs;
tol_iqs      = npar.tol_iqs;

npar.theta_old=[];

% save values at beginning of macro time step: they are needed in the IQS iteration
X_beg=X;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
shape_beg=u_shape(1:npar.n,:);
shape_end=shape_beg(:,end);

tn = tn - dt;
time_beg = tn(end);
time_end = time_beg+dt;

% storage for temp SDIRK quantities
K=zeros(length(shape_beg(:,1)),rk.s);

NFId_old = assemble_mass(dat.nusigf_d,time_beg) / npar.keff;

for iter = 1: max_iter_iqs
    
    % solve for amplitude function
    if strcmpi(npar.prke_solve,'matlab')
        [X,w,t,y]    = solve_prke_ode(X_beg(:,end),dt_macro,time_end,shape_beg(:,end),shape_end);
    else
        [X,dpdt,t,y] = solve_prke_iqs(X_beg(:,end),dt_macro,time_end,shape_beg(:,end),shape_end,npar.n_micro,npar.freq_react);
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
    
    % extract details from piece-wise polynomial by breaking it apart
    [breaks,coefs,l,k,d] = unmkpp(pp);
    % make the polynomial that describes the derivative
    pp_der = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
    
    p_beg    = ppval(pp    ,time_beg);
    dpdt_beg = ppval(pp_der,time_beg);
   
    if strcmp(npar.an_interp_type,'lagrange')
        a = zeros(order+1,2);
    elseif strcmp(npar.an_interp_type,'hermite')
        i2 = zeros(npar.n,2); i3=i2;
        D    = assemble_stiffness(dat.cdiff   ,time_beg);
        A    = assemble_mass(     dat.siga    ,time_beg);
        Aiqs = assemble_mass(     dat.inv_vel ,time_beg) * dpdt_beg/p_beg;
        NFIp = assemble_mass(     dat.nusigf_p,time_beg) / npar.keff;
        IV   = assemble_mass(     dat.inv_vel ,time_beg);
        TR  = NFIp-(D+A+Aiqs);
        rhs = lambda.*C_old/p_beg;
        if npar.set_bc_last
            [TR]=apply_BC_mat_only(TR,npar.add_zero_on_diagonal);
        else
            rhs=apply_BC_vec_only(rhs);
        end
        dshape_beg = IV\(TR*shape_beg(:,end) + rhs);
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:rk.s
        % compute stage time
        ti = time_beg + rk.c(i)*dt;
        t2=ti;
        t1=time_beg;
        t_all = [tn t2];
        dti=rk.c(i)*dt;
        
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
        NFId_new = assemble_mass(dat.nusigf_d ,ti)    / npar.keff ;
        S    = assemble_source(   dat.source_phi,ti)/pi;
        
        % flux-flux matrix
        TR=NFIp-(D+A+Aiqs);
        % build rhs
        zi =  lambda *( C_old*exp(-lambda*dti) )/pi;
        
        if strcmp(npar.an_interp_type,'lagrange')
            for j=1:(order+1)
                Aj = @(t) LagrangeInterp_section(t,t_all,j).*exp(-lambda*(t2-t)).*ppval(pp,t);
                a(j,1) = integral(@(t)Aj(t).*(t2-t)/dti,t1,t2);
                a(j,2) = integral(@(t)Aj(t).*(t-t1)/dti,t1,t2);
            end

            TR = TR + lambda * ( a(end,1)*NFId_old + a(end,2)*NFId_new )/pi ;

            for n=1:order
                zi = zi + lambda *( a(n,1)*NFId_old*shape_beg(:,n) + a(n,2)*NFId_new*shape_beg(:,n) )/pi;
            end

        elseif strcmp(npar.an_interp_type,'hermite')
            I1 = @(t) (t-t1).^2.*ppval(pp,t); %t.^2-2*t1*t+t1^2;
            I2 = @(t) (t-t1).*ppval(pp,t);
            I3 = @(t) ppval(pp,t);
            i1(1)   = integral(@(t)I1(t).*(t2-t)/dti.*exp(-lambda*(t2-t)),t1,t2);
            i2(:,1) = integral(@(t)I2(t).*(t2-t)/dti.*exp(-lambda*(t2-t)),t1,t2)*dshape_beg;
            i3(:,1) = integral(@(t)I3(t).*(t2-t)/dti.*exp(-lambda*(t2-t)),t1,t2)*shape_beg(:,end);
            i1(2)   = integral(@(t)I1(t).*(t-t1)/dti.*exp(-lambda*(t2-t)),t1,t2);
            i2(:,2) = integral(@(t)I2(t).*(t-t1)/dti.*exp(-lambda*(t2-t)),t1,t2)*dshape_beg;
            i3(:,2) = integral(@(t)I3(t).*(t-t1)/dti.*exp(-lambda*(t2-t)),t1,t2)*shape_beg(:,end);
            num = -I2(t2)*dshape_beg - I3(t2)*shape_beg(:,end); % (-dun*t2-un+dun*t1);
            deno = I1(t2); % (t2^2-2*t1*t2+t1^2);

            TR = TR + lambda * ( i1(1)*NFId_old + i1(2)*NFId_new )/deno;
            zi = zi + lambda * ( NFId_old*(i1(1)*num/deno+i2(:,1)+i3(:,1)) + NFId_new*(i1(2)*num/deno+i2(:,2)+i3(:,2)) )/pi;

        end
        zi = zi + S;
        
        if strcmp(npar.method,'BDF')
            % build system matrix
            A = IV - dt_macro*bdf.b(order)*TR;
            % build rhs
            rhs = 0;
            for j=1:order
                rhs = rhs + bdf.a(order,j)*shape_beg(:,j);
            end
            rhs = IV*rhs + dt_macro*bdf.b(order)*zi;
        else    
            % build system matrix
            A = IV-rk.a(i,i)*dt*TR;
            % build rhs
            rhs = IV*shape_beg(:,end) + rk.a(i,i)*dt *zi;
            for j=1:rk.s-1
                rhs = rhs + rk.a(i,j)*dt*K(:,j);
            end
        end
        
        % apply BC
        if npar.set_bc_last
            [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
        else
            rhs=apply_BC_vec_only(rhs);
        end
        % solve: M(unew-uold)/dt=TR.unew
        shape_end = A\rhs;
        % store for temp SDIRK quantities
        K(:,i)=TR*shape_end + zi;
    end
    
    dat.ode.f_end=IV\K(1:npar.n,rk.s);
    
    % update precursors
    if npar.hermite_prec_update
        mat=[ t1^3 t1^2 t1 1; t2^3 t2^2 t2 1;  3*t1^2 2*t1 1 0; 3*t2^2 2*t2 1 0];
        rhs = [ shape_beg(:,end)' ;shape_end' ; dat.ode.f_beg'; dat.ode.f_end'];
        % contains the w coefficients such that :
        %  Phi(t) = w(1) t^3 + w(2) t^2 + w(3) t + w(4)
        w = mat\rhs; w=w'; nw=size(w,2);
        b = zeros(nw,2);
        for k=1:nw
            b(k,1) = integral(@(t) (t2-t)/dt .*t.^(nw-k) .*exp(-lambda*(t2-t)).*ppval(pp,t),t1,t2);
            b(k,2) = integral(@(t) (t-t1)/dt .*t.^(nw-k) .*exp(-lambda*(t2-t)).*ppval(pp,t),t1,t2);
        end
        tmp = (w*b);
        C_new =  C_old*exp(-lambda*dt) + NFId_old*tmp(:,1) + NFId_new*tmp(:,2);
    elseif strcmp(npar.an_interp_type,'lagrange')
        C_new =  C_old*exp(-lambda*dt) + ( a(end,1)*NFId_old + a(end,2)*NFId_new )*shape_end ;
        for n=1:order
            C_new = C_new + ( a(n,1)*NFId_old + a(n,2)*NFId_new )*shape_beg(:,n);
        end
    elseif strcmp(npar.an_interp_type,'hermite')
        C_new = C_old*exp(-lambda*dt) + ( (i1(1)*NFId_old + i1(2)*NFId_new )/deno )*shape_end  + ( NFId_old*(i1(1)*num/deno+i2(:,1)+i3(:,1)) + NFId_new*(i1(2)*num/deno+i2(:,2)+i3(:,2)) );
    else
        error('Interpolation type must be lagrange or hermite')
    end
    
    % re-package as single solution vector
    u_shape = [ shape_end ; C_new];
    
    % check for tolerance
    err = abs( (npar.phi_adj)'*IV*shape_end/npar.K0  - 1);
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
%     warning('IQS did not converge in %s',mfilename);
end


% renormalize anyway
u_shape = u_shape / ( ((npar.phi_adj)'*npar.IV*shape_end) / npar.K0 );