function [ u ] = solve_TD_diffusion_JFNK(u,dt,tn)

global io dat npar

lambda = dat.lambda;

% split initial solution vector into flux and precursors:
Phi_old = u(1:npar.n,:);
Phi_new = u(1:npar.n,end);
C_old   = u(npar.n+1:end,:);

order = length(tn) - 1;
time_end = tn(end);
time_beg = time_end-dt;

% assemble matrices
D    = assemble_stiffness(dat.cdiff   ,time_end);
A    = assemble_mass(     dat.siga    ,time_end);
NFIp = assemble_mass(     dat.nusigf_p,time_end) / npar.keff;
IV   = assemble_mass(     dat.inv_vel ,time_end);
NFId_new = assemble_mass(dat.nusigf_d ,time_end) / npar.keff;
NFId_old = assemble_mass(dat.nusigf_d ,time_beg) / npar.keff;

TR=NFIp -(D+A);

% Ci = C_old; 
% deno = (1+lambda*dt);
% Ci = Ci/deno;
% 
% TR = TR + lambda*dt/deno*NFId_new;
% zi = lambda*Ci;

zi =  lambda*(C_old*exp(-lambda*dt));

for j=1:(order+1)
    Aj = @(t) LagrangeInterp_section(t,tn,j).*exp(-lambda*(time_end-t));
    a(j,1) = integral(@(t)Aj(t).*(time_end-t)/dt,time_beg,time_end);
    a(j,2) = integral(@(t)Aj(t).*(t-time_beg)/dt,time_beg,time_end);
end

TR = TR + lambda * ( a(end,1)*NFId_old + a(end,2)*NFId_new ) ;

for n=1:order
    zi = zi + lambda *( a(n,1)*NFId_old*Phi_old(:,n) + a(n,2)*NFId_new*Phi_old(:,n) );
end

A = IV - dt*TR;
rhs = IV*Phi_old(:,end) + dt*zi;

[A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);

Re_fun = @(phi) A*phi - rhs;

eps = 1e-6;
Jac_fun = @(phi,y) (Re_fun(phi+eps*y) - Re_fun(phi)) / eps; 

% Compute residual
Re = Re_fun(Phi_new);
residual = norm(Re,2);
if io.console_print
    fprintf(' Nonlinear|R| %g \n',residual);
end

% Solve linear system
dphi = gmres(@(y) Jac_fun(Phi_new,y), -Re,[],1e-4,npar.n);
Phi_new = Phi_new + dphi;

% C_new = Ci +  dt / deno * NFId_new*Phi_new;
C_new =  C_old*exp(-lambda*dt) + ( a(end,1)*NFId_old + a(end,2)*NFId_new )*Phi_new ;
for n=1:order
    C_new = C_new + ( a(n,1)*NFId_old + a(n,2)*NFId_new )*Phi_old(:,n);
end
% re-package as single solution vector
u = [Phi_new ; C_new];
