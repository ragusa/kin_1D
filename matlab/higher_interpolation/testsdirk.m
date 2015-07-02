clear all; clc;
% du/dt = u + q

g=0.43586652150845899941601945119356;
A=[g 0 0; ...
    ((1-g)/2) g 0;...
    (-(6*g^2-16*g+1)/4) ((6*g^2-20*g+5)/4) g];
c=sum(A'); b=A(end,:);

% if mms u is quadratic, u=1+t+t^2, thus q=t-t^2
u = @(t)(1+t+t^2);
q = @(t)(t*(1-t));
% if mms u is linear, use this:
% u = @(t)(1+t);
% q = @(t)(-t);

yn=u(0); 
dt=1.;
time0=0;

% stage 1
t1 = time0 + c(1)*dt;
deno = 1 - dt*A(1,1);
Y1 = ( yn + dt*A(1,1)*q(t1) )/deno;
% stage 2
t2 = time0 + c(2)*dt;
deno = 1 - dt*A(2,2);
Y2 = ( yn + dt*A(2,1)*(Y1+q(t1)) + dt*A(2,2)*q(t2) )/deno;
% stage 3
t3 = time0 + c(3)*dt;
deno = 1 - dt*A(3,3);
Y3 = ( yn + dt*A(3,1)*(Y1+q(t1)) + dt*A(3,2)*(Y2+q(t2)) + dt*A(3,3)*q(t3) )/deno;

% printout
fprintf('end time step num value %15.10g\n',Y3);
fprintf('end time step exa value %15.10g\n',u(dt));
fprintf('difference %15.10g\n',Y3-u(dt));
