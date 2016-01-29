function [ u ] = solve_TD_diffusion_an_prec(u,dt,time_end)

global dat npar

% split initial solution vector into flux and precursors:
Phi_old = u(1:npar.n);
C_old   = u(npar.n+1:end);

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

% beginning of the time interval
time_beg = time_end - dt;
ti = time_beg;
% storage for temp SDIRK quantities
K=zeros(length(Phi_old),rk.s);

Phi = cell(rk.s+1,1);
Phi{1} = Phi_old;
NFId = cell(rk.s+1,1);
NFId{1} = assemble_mass(dat.nusigf_d,time_beg) / npar.keff;

for i=1:rk.s
    Phi{i+1} = zeros(size(Phi{i}));
    
    % compute stage time
    ti_old = ti;
    ti = time_beg + rk.c(i)*dt;
    
    % build system matrix
    D    = assemble_stiffness(dat.cdiff   ,ti);
    A    = assemble_mass(     dat.siga    ,ti);
    NFIp = assemble_mass(     dat.nusigf_p,ti) / npar.keff;
    IV   = assemble_mass(     dat.inv_vel ,ti);
    NFId{i+1} = assemble_mass( dat.nusigf_d,ti) / npar.keff;
    
    % flux-flux matrix
    TR=NFIp-(D+A);
    
    % integrals for the analytical expressions for the precursors
    t2=ti;
    t1=ti_old;
    dti=t2-t1;
    
    A1= @(t)( ((t2-t)/dti).^2      .*exp(-lambda*(t2-t)) );
    A2= @(t)( (t-t1).*(t2-t)/dti^2 .*exp(-lambda*(t2-t)) );
    A3= @(t)( ((t-t1)/dti).^2      .*exp(-lambda*(t2-t)) );
    
    a1(i)= integral(@(t)A1(t),t1,t2);
    a2(i)= integral(@(t)A2(t),t1,t2);
    a3(i)= integral(@(t)A3(t),t1,t2);
    
    % build transient matrix
    TR = TR + lambda * ( a2(i)*NFId{i} + a3(i)*NFId{i+1} ) ;
    
    % build system matrix
    A = IV-rk.a(i,i)*dt*TR;
    
    % build rhs
    zi =  lambda *( C_old*exp(-lambda*dt*rk.c(i)) );
    for k=1:i
        zi = zi + lambda *( ( a1(k)*NFId{k} + a2(k)*NFId{k+1} )*Phi{k} + ( a2(k)*NFId{k} + a3(k)*NFId{k+1} )*Phi{k+1} );
    end
    rhs = IV*Phi_old + rk.a(i,i)*dt *zi;
    for j=1:rk.s-1
        rhs = rhs + rk.a(i,j)*dt*K(:,j);
    end
    % apply BC
    if npar.set_bc_last
        [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
    else
        rhs=apply_BC_vec_only(rhs);
    end
    % solve: M(unew-uold)/dt=TR.unew
    Phi{i+1} = A\rhs;
    % store for temp SDIRK quantities
    K(:,i)=TR*Phi{i+1} + zi;
end
Phi_new = Phi{end};

% update precursors
if npar.rk.s==1;
    C_new =  C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*Phi_old + ( a2*NFId_old + a3*NFId_new )*Phi_new ;
else
    if npar.iqs_prke_interpolation_method>=3
        t1 = time_beg; t2 = time_end;
        NFId_old = NFId{1}; NFId_new = NFId{end};
        dat.ode.f_end=K(1:npar.n,rk.s);
        mat=[ t1^3 t1^2 t1 1; t2^3 t2^2 t2 1;  3*t1^2 2*t1 1 0; 3*t2^2 2*t2 1 0];
        rhs = [ Phi_old' ;Phi_new' ; dat.ode.f_beg'; dat.ode.f_end'];
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
    else
        C_new =  C_old*exp(-lambda*dt) + ( a1(end)*NFId_old + a2(end)*NFId_new )*Phi_old + ( a2(end)*NFId_old + a3(end)*NFId_new )*Phi_new ;        
    end
end
% re-package as single solution vector
u = [Phi_new;C_new];


end

