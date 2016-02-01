function [out] = SDIRK_an(yn,tn,A,a_rk,c_rk,dt)

[var,order] = size(yn);
if length(tn) ~= order
    error('Number of solutions must match number of times')
end
a = zeros(order+1,1);
stages = length(c_rk);
k=zeros(1,stages);
Y=zeros(var,stages);
tt = tn(end)+c_rk*dt;
t_all = [tn 0];

for i=1:stages
    t_all(end) = tt(i);
    Ai  = A(tt(i));
    dti = c_rk(i)*dt;
    
    for j=1:(order+1)
        a(j) = integral(@(t) LagrangeInterp_section(t,t_all,j).*exp(Ai(2,2)*(tt(i)-t)), tn(end),tt(i));
    end
    
    Ki = Ai(1,1) + Ai(1,2)*Ai(2,1)*a(end);
    
    zi = Ai(1,2)*yn(2,end)*exp(Ai(2,2)*(tt(i)-tn(end)));
    for n=1:order
        zi = zi + Ai(1,2)*Ai(2,1)*yn(1,n)*a(n);
    end
    
    M = 1-a_rk(i,i)*dt*Ki;
    rhs = yn(1,end) + a_rk(i,i)*dt*zi;
    for j=1:i-1,
        rhs = rhs + dt*a_rk(i,j)*k(:,j);
    end
    Y(1,i) = M\rhs;
    k(i) = Ki*Y(1,i) + zi;
end
Y(2,end) = yn(2,end)*exp(Ai(2,2)*(tt(i)-tn(end))) + Ai(2,1)*a(end)*Y(1,end);
for n=1:order
    Y(2,end) = Y(2,end) + Ai(2,1)*yn(1,n)*a(n);
end
out=Y(:,end);