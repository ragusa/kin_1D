function [ u ] = solve_TD_diffusion_an_prec(u,dt,time_end)

global dat npar

% split initial solution vector into flux and precursors:
Phi_old = u(1:npar.n);
C_old   = u(npar.n+1:end);

% shortcut
lambda = dat.lambda;

% (Phi^new - Phi^old)/dt = {[Fisp - Leak - Abso] Phi }^new + lambda C^new
%
% C^new = C^old exp(-lambda dt)   + exp(-lambda*t^end) ...
%   * int_{t^beg}^{t^end} exp(lambda.t) Fisd(t) Phi(t) dt
% Note that:
% exp(-lambda*t^end) * int_{t^beg}^{t^end} exp(lambda.t) Fisd(t) Phi(t) dt
% = int_{t^beg}^{t^end} exp(-lambda.(t^end-t) Fisd(t) Phi(t) dt
% = {Fisd.Phi}^old a1 + {Fisd^old.Phi^new + Fisd^new.Phi^old} a2 + {Fisd.Phi}^new a3 
% with
% a1 = int_{t^beg}^{t^end}exp(-lambda.(t^end-t) [(t2-t)/dt]^2      dt
% a2 = int_{t^beg}^{t^end}exp(-lambda.(t^end-t) (t-t1)(t2-t)/dt^2  dt
% a3 = int_{t^beg}^{t^end}exp(-lambda.(t^end-t) [(t-t1)/dt]^2      dt

TR  = sparse(npar.n,npar.n,npar.nnz);

D    = assemble_stiffness(dat.cdiff   ,time_end);
A    = assemble_mass(     dat.siga    ,time_end);
NFIp = assemble_mass(     dat.nusigf_p,time_end) / npar.keff;
IV   = assemble_mass(     dat.inv_vel ,time_end);

NFId_old = assemble_mass(dat.nusigf_d,time_end-dt) / npar.keff;
NFId_new = assemble_mass(dat.nusigf_d,time_end)    / npar.keff;

% flux-flux matrix
TR=NFIp-(D+A);

% integrals for the analytical expressions for the precursors
t2=time_end;
t1=t2-dt;

A1= @(t)( ((t2-t)/dt).^2      .*exp(-lambda*(t2-t)) );
A2= @(t)( (t-t1).*(t2-t)/dt^2 .*exp(-lambda*(t2-t)) );
A3= @(t)( ((t-t1)/dt).^2      .*exp(-lambda*(t2-t)) );

a1= quad(@(t)A1(t),t1,t2);
a2= quad(@(t)A2(t),t1,t2);
a3= quad(@(t)A3(t),t1,t2);


% build transient matrix
TR = TR + lambda * ( a2*NFId_old + a3*NFId_new ) ;

% build rhs from backward Euler time discretization
rhs = IV*Phi_old + dt * lambda *( C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*Phi_old );

% build system matrix
A = IV-dt*TR;

% apply BC
if npar.set_bc_last
    [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
else
    rhs=apply_BC_vec_only(rhs);
end

% solve
Phi_new = A\rhs;

% update precursors
C_new =  C_old*exp(-lambda*dt) + ( a1*NFId_old + a2*NFId_new )*Phi_old + ( a2*NFId_old + a3*NFId_new )*Phi_new ;

% re-package as single solution vector
u = [Phi_new;C_new];



end

