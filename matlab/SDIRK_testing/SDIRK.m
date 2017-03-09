function [out] = SDIRK(yn,tn,A,a,c,dt)

stages = length(c);
var = length(yn);
k=zeros(var,stages);
Y=zeros(var,stages);
t = tn+c*dt;

for i=1:stages
    M = eye(var)-a(i,i)*dt*A(t(i));
    rhs = yn;
    for j=1:i-1,
        rhs = rhs + dt*a(i,j)*k(:,j);
    end
    Y(:,i) = M\rhs;
    X = Y(:,i);
    k(:,i) = A(t(i))*Y(:,i);
end
out=Y(:,end);
