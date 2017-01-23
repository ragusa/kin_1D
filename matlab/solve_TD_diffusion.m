function [ u ] = solve_TD_diffusion(u0,dt,tn)

global npar

% % TR = assemble_transient_operator(time_end);
% % M  = assemble_time_dependent_operator(time_end);
% % 
% % % build rhs from backward Euler time discretization
% % rhs = M*u;
% % % build system matrix
% % A = M-dt*TR;
% % if npar.set_bc_last
% %     [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
% % else
% %     rhs=apply_BC_vec_only(rhs);
% % end
% % % solve: M(unew-uold)/dt=TR.unew
% % u = A\rhs;

% SDIRK solve
% Yi = yn + h sum_j aij f(tj,Yj)
% gives (f=Qy)
%  [I-aii.h.Q(ti)]Yi = yn + h sum(j=1..i-1) { aij Q(tj)Yj }
rk = npar.rk;
if strcmp(npar.method,'BDF')
    bdf = npar.bdf;
    order = length(tn);
end

% beginning of the time interval
time_beg = tn(end) - dt;
% storage for temp SDIRK quantities
K=zeros(length(u0),rk.s-1);

for i=1:rk.s
    % compute stage time
    ti = time_beg + rk.c(i)*dt;
    TR = assemble_transient_operator(ti);
    M  = assemble_time_dependent_operator(ti);
    S = assemble_source_operator(ti);
    
    if strcmp(npar.method,'BDF')
        % build system matrix
        A = M - dt*bdf.b(order)*TR;
        % build rhs
        rhs = 0;
        for j=1:order
            rhs = rhs + bdf.a(order,j)*u0(:,j);
        end
        rhs = M*rhs + dt*bdf.b(order)*S;

    else
        % build system matrix
        A = M-rk.a(i,i)*dt*TR;
        % build rhs
        rhs = M*u0 + rk.a(i,i)*dt*S;
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
    u = A\rhs;
    % store for temp SDIRK quantities
    if i<rk.s
        K(:,i)=TR*u + S;
    end
end

end

