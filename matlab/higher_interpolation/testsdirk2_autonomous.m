function testsdirk2_autonomous
clc
% du/dt = u + q
% 
% SDIRK33
g=0.43586652150845899941601945119356;
A=[g 0 0; ...
    ((1-g)/2) g 0;...
    (-(6*g^2-16*g+1)/4) ((6*g^2-20*g+5)/4) g];
% SDIRK54
% sdirk54 constants:
A=[ 1./4., 0., 0., 0., 0.;...
    1./2., 1./4., 0., 0., 0.;...
    17./50., -1./25., 1./4., 0., 0.;...
    371./1360., -137./2720., 15./544., 1./4., 0.;...
    25./24., -49./48., 125./16., -85./12., 1./4.];
% common to SDIRK
c=sum(A'); b=A(end,:);

% if mms u is quadratic, u=1+t+t^2, thus q=t-t^2
u = @(t)(1+t+t.^2);
q = @(t)(t.*(1-t));
% if mms u is linear, use this:
% u = @(t)(1+t);
% q = @(t)(-t);

yn=[u(0);0];
dt=1e-1;
time0=0;

% number of stages
n_stages = length(c);
Y=zeros(2,n_stages); F=Y;
% Yi = yn + dt sum_j { A_ij f(tj, Yj) }
for i=1:n_stages
    % stage time
    ts(i) = time0 + c(i)*dt;
    % denominator
    aux = yn;
    for j=1:i-1
        aux = aux + dt*A(i,j)*F(:,j);
    end
    [Y(:,i),F(:,i)] = my_solve(Y(:,i),q,aux,dt*A(i,i));
end

Y
% printout
fprintf('end time step num value %15.10g\n',Y(1,end));
fprintf('end time step exa value %15.10g\n',u(dt));
fprintf('difference %15.10g\n',Y(1,end)-u(dt));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,F]=my_solve(X,q,aux,dti)

% aux = yn + sum_{j<i} dt*aij*Fj
% dti = dt*aii

for iter=1:100
    Xnew = aux + dti * ss_residual(X,q);
    if abs(Xnew-X)<1e-12
        X=Xnew;
        break
    end
    if iter==100
        [Xnew X]
        error('not converged');
    end
    X=Xnew;
end

F=ss_residual(X,q);

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F=ss_residual(X,q)

F=zeros(2,1);
F(1) = X(1) + q(X(2));
F(2) = 1;

return
end
