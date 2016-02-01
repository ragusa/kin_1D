function [out] = SDIRK_an(yn,tn,A,a,c,dt)

stages = length(c);
var = length(yn);
k=zeros(1,stages);
Y=zeros(var,stages);
tt = tn+c*dt;

for i=1:stages
    Ai  = A(tt(i));
    dti = c(i)*dt;
    a1 = integral(@(t) (tt(i)-t)/dti .*exp(Ai(2,2)*(tt(i)-t)), tn,tt(i));
    a2 = integral(@(t) (t-tn)/dti    .*exp(Ai(2,2)*(tt(i)-t)), tn,tt(i));
    Ki = Ai(1,1) + Ai(1,2)*Ai(2,1)*a2;
    zi = Ai(1,2)*( yn(2)*exp(Ai(2,2)*(tt(i)-tn)) + Ai(2,1)*a1*yn(1) );
    M = 1-a(i,i)*dt*Ki;
    rhs = yn(1) + a(i,i)*dt*zi;
    for j=1:i-1,
        rhs = rhs + dt*a(i,j)*k(:,j);
    end
    Y(1,i) = M\rhs;
    k(i) = Ki*Y(1,i) + zi;
end
Y(2,end) = yn(2)*exp(Ai(2,2)*(tt(i)-tn)) + Ai(2,1)*(a1*yn(1) + a2*Y(1,end));
out=Y(:,end);