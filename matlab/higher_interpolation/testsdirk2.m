clear all; clc;
% du/dt = u + q

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
u = @(t)(t.^2);
q = @(t)(t.*(2-t));
% % if mms u is qintic, u=t^5, thus q=5*t^4-t^5
% u = @(t)(t.^5);
% q = @(t)(t.^4.*(5-t));
% if mms u is linear, use this:
% u = @(t)(1+t);
% q = @(t)(-t);

yn=u(0);
dt=1;
time0=0;

% number of stages
n_stages = length(c);
Y=zeros(n_stages,1); ts=Y;
% Yi = yn + dt sum_j { A_ij f(tj, Yj) }
for i=1:n_stages
    % stage time
    ts(i) = time0 + c(i)*dt;
    % denominator
    deno = 1 - dt*A(i,i);
    % numerator
    nume = yn + dt*A(i,i)*q(ts(i)) ;
    for j=1:i-1
        nume = nume + dt*A(i,j)*(Y(j)+q(ts(j)));
    end
    Y(i) = nume/deno;
end
Y
% printout
fprintf('end time step num value %15.10g\n',Y(end));
fprintf('end time step exa value %15.10g\n',u(dt));
fprintf('relative difference %15.10g\n',Y(end)/u(dt)-1);

% other approach
M=eye(n_stages) - dt*A;
rhs = yn*ones(n_stages,1) + dt*A*q(ts);
YY=M\rhs

iM=inv(M);
r=iM(end,:);
dot(r, yn+dt*A*q(ts))


k=zeros(n_stages,1); ts=Y;
% ki = f( ti, yn + dt sum_j { A_ij kj) }
for i=1:n_stages
    % stage time
    ts(i) = time0 + c(i)*dt;
    % denominator
    deno = 1 - dt*A(i,i);
    % numerator
    nume = yn + q(ts(i));
    for j=1:i-1
        nume = nume + dt*A(i,j)*k(j);
    end
    k(i) = nume/deno;
end
ynew = yn + dt* dot(b,k);
fprintf('end time step num value %15.10g\n',ynew);
fprintf('end time step exa value %15.10g\n',u(dt));
fprintf('relative difference %15.10g\n',ynew/u(dt)-1);

fprintf('\ndifference between Y and k %15.10g\n',ynew-Y(end));
 
[k, u(ts)+q(ts)]