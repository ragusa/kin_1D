function [ varargout ] = solve_TD_diffusion_an_prec(u,dt,tn)

global dat npar

% split initial solution vector into flux and precursors:
Phi_old = u(1:npar.n,:);
C_old   = u(npar.n+1:end,end);
tn = tn-dt;

% shortcut
lambda = dat.lambda;
% %
% % % (Phi^new - Phi^old)/dt = {[Fisp - Leak - Abso] Phi }^new + lambda C^new
% % %
% % % C^new = C^old exp(-lambda dt)   + exp(-lambda*t^end) ...
% % %   * int_{t^beg}^{t^end} exp(lambda.t) Fisd(t) Phi(t) dt
% % % Note that:
% % % exp(-lambda*t^end) * int_{t^beg}^{t^end} exp(lambda.t) Fisd(t) Phi(t) dt
% % % = int_{t^beg}^{t^end} exp(-lambda.(t^end-t) Fisd(t) Phi(t) dt
% % % = {Fisd.Phi}^old a1 + {Fisd^old.Phi^new + Fisd^new.Phi^old} a2 + {Fisd.Phi}^new a3
% % % with
% % % a1 = int_{t^beg}^{t^end}exp(-lambda.(t^end-t) [(t2-t)/dt]^2      dt
% % % a2 = int_{t^beg}^{t^end}exp(-lambda.(t^end-t) (t-t1)(t2-t)/dt^2  dt
% % % a3 = int_{t^beg}^{t^end}exp(-lambda.(t^end-t) [(t-t1)/dt]^2      dt
% %
% % TR  = sparse(npar.n,npar.n,npar.nnz);
% %
% % D    = assemble_stiffness(dat.cdiff   ,time_end);
% % A    = assemble_mass(     dat.siga    ,time_end);
% % NFIp = assemble_mass(     dat.nusigf_p,time_end) / npar.keff;
% % IV   = assemble_mass(     dat.inv_vel ,time_end);
% %
% % NFId_old = assemble_mass(dat.nusigf_d,time_end-dt) / npar.keff;
% % NFId_new = assemble_mass(dat.nusigf_d,time_end)    / npar.keff;
% %
% % % flux-flux matrix
% % TR=NFIp-(D+A);
% %
% % % integrals for the analytical expressions for the precursors
% % t2=time_end;
% % t1=t2-dt;
% %
% % A1= @(t)( ((t2-t)/dt).^2      .*exp(-lambda*(t2-t)) );
% % A2= @(t)( (t-t1).*(t2-t)/dt^2 .*exp(-lambda*(t2-t)) );
% % A3= @(t)( ((t-t1)/dt).^2      .*exp(-lambda*(t2-t)) );
% %
% % a1= integral(@(t)A1(t),t1,t2);
% % a2= integral(@(t)A2(t),t1,t2);
% % a3= integral(@(t)A3(t),t1,t2);
% %
% %
% % % build transient matrix
% % TR = TR + lambda * ( a2*NFId_old + a3*NFId_new ) ;
% %
% % % build rhs from backward Euler time discretization
% % rhs = IV*Phi_old + dt * lambda *( C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*Phi_old );
% %
% % % build system matrix
% % A = IV-dt*TR;
% %
% % % apply BC
% % if npar.set_bc_last
% %     [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
% % else
% %     rhs=apply_BC_vec_only(rhs);
% % end
% %
% % % solve
% % Phi_new = A\rhs;
% %
% % % update precursors
% % C_new =  C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*Phi_old + ( a2*NFId_old + a3*NFId_new )*Phi_new ;
% %
% % % re-package as single solution vector
% % u = [Phi_new;C_new];
% %
% % end




% SDIRK solve
% Yi = yn + h sum_j aij f(tj,Yj)
% gives (f=Qy+z)
%  [I-aii.h.Q(ti)]Yi = yn + h.aii zi + h sum(j=1..i-1) ( aij {Q(tj)Yj + zj} )
rk=npar.rk;
if strcmp(npar.method,'BDF')
    bdf = npar.bdf;
end

% beginning of the time interval
time_beg = tn(end);
order = length(tn);

if strcmp(npar.an_interp_type,'lagrange')
    a = zeros(order+1,2);
elseif strcmp(npar.an_interp_type,'hermite')
    i2 = zeros(npar.n,2); i3=i2;
    D    = assemble_stiffness(dat.cdiff   ,time_beg);
    A    = assemble_mass(     dat.siga    ,time_beg);
    NFIp = assemble_mass(     dat.nusigf_p,time_beg) / npar.keff;
    IV   = assemble_mass(     dat.inv_vel ,time_beg);
    TR  = (NFIp-(D+A));
    rhs = lambda.*C_old;
    if npar.set_bc_last
        [TR]=apply_BC_mat_only(TR,npar.add_zero_on_diagonal);
    else
        rhs=apply_BC_vec_only(rhs);
    end
    
    dPhi_old = IV\(TR*Phi_old(:,end) + rhs);
end

% storage for temp SDIRK quantities
K=zeros(length(Phi_old(:,1)),rk.s);

NFId_old = assemble_mass(dat.nusigf_d,time_beg) / npar.keff;
auxPhi=zeros(length(Phi_old(:,end)),rk.s);

for i=1:rk.s
    
    % compute stage time
    ti = time_beg + rk.c(i)*dt;
    
    % build system matrix
    D    = assemble_stiffness(dat.cdiff   ,ti);
    A    = assemble_mass(     dat.siga    ,ti);
    NFIp = assemble_mass(     dat.nusigf_p,ti) / npar.keff;
    IV   = assemble_mass(     dat.inv_vel ,ti);
    NFId_new = assemble_mass( dat.nusigf_d,ti) / npar.keff;
    S    = assemble_source(   dat.source_phi,ti);
    
    % integrals for the analytical expressions for the precursors
    t2=ti;
    t1=time_beg;
    t_all = [tn t2];
    dti=rk.c(i)*dt;
    
    % flux-flux matrix
    TR=NFIp-(D+A);        
    % build rhs
    zi =  lambda*(C_old*exp(-lambda*dti)) + S;
    
    if strcmp(npar.an_interp_type,'lagrange')
        for j=1:(order+1)
            Aj = @(t) LagrangeInterp_section(t,t_all,j).*exp(-lambda*(t2-t));
            a(j,1) = integral(@(t)Aj(t).*(t2-t)/dti,t1,t2);
            a(j,2) = integral(@(t)Aj(t).*(t-t1)/dti,t1,t2);
        end
    
        TR = TR + lambda * ( a(end,1)*NFId_old + a(end,2)*NFId_new ) ;
    
        for n=1:order
            zi = zi + lambda *( a(n,1)*NFId_old*Phi_old(:,n) + a(n,2)*NFId_new*Phi_old(:,n) );
        end
        
    elseif strcmp(npar.an_interp_type,'hermite')
        I1 = @(t) (t-t1).^2; %t.^2-2*t1*t+t1^2;
        I2 = @(t) (t-t1);
        I3 = @(t) 1;
        i1(1)   = integral(@(t)I1(t).*(t2-t)/dti.*exp(-lambda*(t2-t)),t1,t2);
        i2(:,1) = integral(@(t)I2(t).*(t2-t)/dti.*exp(-lambda*(t2-t)),t1,t2)*dPhi_old;
        i3(:,1) = integral(@(t)I3(t).*(t2-t)/dti.*exp(-lambda*(t2-t)),t1,t2)*Phi_old(:,end);
        i1(2)   = integral(@(t)I1(t).*(t-t1)/dti.*exp(-lambda*(t2-t)),t1,t2);
        i2(:,2) = integral(@(t)I2(t).*(t-t1)/dti.*exp(-lambda*(t2-t)),t1,t2)*dPhi_old;
        i3(:,2) = integral(@(t)I3(t).*(t-t1)/dti.*exp(-lambda*(t2-t)),t1,t2)*Phi_old(:,end);
        num = -I2(t2)*dPhi_old - I3(t2)*Phi_old; % (-dun*t2-un+dun*t1);
        deno = I1(t2); % (t2^2-2*t1*t2+t1^2);
        
        TR = TR + lambda * ( i1(1)*NFId_old + i1(2)*NFId_new )/deno;
        zi = zi + lambda * ( NFId_old*(i1(1)*num/deno+i2(:,1)+i3(:,1)) + NFId_new*(i1(2)*num/deno+i2(:,2)+i3(:,2)) );
        
    end
    
    if strcmp(npar.method,'BDF')
        % build system matrix
        A = IV - dt*bdf.b(order)*TR;
        % build rhs
        rhs = 0;
        for j=1:order
            rhs = rhs + bdf.a(order,j)*Phi_old(:,j);
        end
        rhs = dt*bdf.b(order)*zi + IV*rhs;        
    else
        % build system matrix
        A = IV-rk.a(i,i)*dt*TR;
        % build rhs
        rhs = IV*Phi_old(:,end) + rk.a(i,i)*dt *zi;
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
    Phi_new = A\rhs;
    % save intermediate fluxes
    auxPhi(:,i)=Phi_new;
    % store for temp SDIRK quantities
    K(:,i)=TR*Phi_new + zi;
end

dat.ode.f_end=IV\K(1:npar.n,rk.s);

% update precursors
if npar.hermite_prec_update
    mat=[ t1^3 t1^2 t1 1; t2^3 t2^2 t2 1;  3*t1^2 2*t1 1 0; 3*t2^2 2*t2 1 0];
    rhs = [ Phi_old(:,end)' ;Phi_new' ; dat.ode.f_beg'; dat.ode.f_end'];
    % contains the w coefficients such that :
    %  Phi(t) = w(1) t^3 + w(2) t^2 + w(3) t + w(4)
    w = mat\rhs; w=w'; nw=size(w,2);
    b = zeros(nw,2);
    for k=1:nw
        b(k,1) = integral(@(t) (t2-t)/dt .*t.^(nw-k) .*exp(-lambda*(t2-t)),t1,t2);
        b(k,2) = integral(@(t) (t-t1)/dt .*t.^(nw-k) .*exp(-lambda*(t2-t)),t1,t2);
    end
    tmp = w*b;
    C_new =  C_old*exp(-lambda*dt) + NFId_old*tmp(:,1) + NFId_new*tmp(:,2);
elseif strcmp(npar.an_interp_type,'lagrange')
    C_new =  C_old*exp(-lambda*dt) + ( a(end,1)*NFId_old + a(end,2)*NFId_new )*Phi_new ;
    for n=1:order
        C_new = C_new + ( a(n,1)*NFId_old + a(n,2)*NFId_new )*Phi_old(:,n);
    end
elseif strcmp(npar.an_interp_type,'hermite')
    C_new = C_old*exp(-lambda*dt) + ( (i1(1)*NFId_old + i1(2)*NFId_new )/deno )*Phi_new  + ( NFId_old*(i1(1)*num/deno+i2(:,1)+i3(:,1)) + NFId_new*(i1(2)*num/deno+i2(:,2)+i3(:,2)) );
else
    error('Interpolation type must be lagrange or hermite')
end
% re-package as single solution vector
u = [Phi_new;C_new];

% output
nOutputs = nargout;
varargout = cell(1,nOutputs);
switch nOutputs
    case(1)
        varargout{1} = u;
    case(2)
        varargout{1} = u;
        varargout{2} = auxPhi;
    otherwise
        error('Wrong number of output arguments in %s',mfilename);
end


end

