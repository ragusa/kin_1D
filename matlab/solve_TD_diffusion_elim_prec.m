function [ varargout ] = solve_TD_diffusion_elim_prec(u,dt,tn)

global dat npar

% split initial solution vector into flux and precursors:
Phi_old = u(1:npar.n,:);
C_old   = u(npar.n+1:end,:);

% shortcut
lambda = dat.lambda;
% %

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
time_beg = tn(end) - dt;
% storage for temp SDIRK quantities (-lambda Ci + Bi Phii)
f=zeros(length(Phi_old(:,end)),rk.s);
fC=f;
% NFId_old = assemble_mass(dat.nusigf_d,time_beg) / npar.keff;
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
    
    % flux-flux matrix
    TR=NFIp-(D+A);
    
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
        zi = lambda*Ci+S;
        rhs = 0;
        for j=1:order
            rhs = rhs + bdf.a(order,j)*Phi_old(:,j);
        end
        rhs = dt*bdf.b(order)*zi + IV*rhs;
    else    
        % expressions for the precursors
        % Ci = [ C_old + h aii FISdi.Phii + sum_{j<i} h aij(-lambda.Cj + FISdj.Phij)]/(1+lambda h aii) 
        % let fci = -lambda.Ci + FISdi.Phii
        Ci = C_old; 
        for j=1:i-1
            Ci = Ci + rk.a(i,j)*dt*fC(:,j);
        end
        deno = (1+lambda*rk.a(i,i)*dt);
        Ci = Ci/deno; % this is Ci without: h.aii.FISdi.Phii/deno

        % build transient matrix
        TR = TR + (lambda * ( rk.a(i,i)*dt / deno )) * NFId_new;     

        % build system matrix
        A = IV-rk.a(i,i)*dt*TR;

        % build rhs
        zi = lambda*Ci + S;
        rhs = IV*Phi_old + rk.a(i,i)*dt *zi;
        for j=1:i-1
            rhs = rhs + rk.a(i,j)*dt*f(:,j);
        end
    end
    
    % apply BC
    if npar.set_bc_last
        [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
    else
        rhs=apply_BC_vec_only(rhs);
    end
    % solve
    Phi_new = A\rhs;
    % finish Ci
    if strcmp(npar.method,'BDF')
        Ci = Ci +  bdf.b(order)*dt / deno * NFId_new*Phi_new;
    else
        Ci = Ci +  rk.a(i,i)*dt / deno * NFId_new*Phi_new;
    end
    % save intermediate fluxes
    auxPhi(:,i)=Phi_new;
    % store for temp SDIRK quantities
    f(:,i)=TR*Phi_new + zi; % caveat: here, TR is: [ NFIp-(D+A) + lambda*aii*dt/deno* NFId_i] so we use zi and not the full Ci
    fC(:,i)=-lambda*Ci + NFId_new*Phi_new;
end

% save for hermite interp
dat.ode.f_end=IV\f(:,rk.s);


% re-package as single solution vector
u = [Phi_new ; Ci];

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

