function [out] = SDIRK_an_herm(yn,tn,A,a_rk,c_rk,dt)

[var,order] = size(yn);
stages = length(c_rk);
k=zeros(1,stages);
Y=zeros(var,stages);
tt = tn+c_rk*dt;
un = yn(1);
dyn = A(tn)*yn; dun = dyn(1);

for i=1:stages
    Ai  = A(tt(i));
    t1 = tn;
    t2 = tt(i);
    dti = c_rk(i)*dt;
    
    I1 = @(t) t.^2-2*t1*t+t1^2;
    I2 = @(t) dun*t+un-dun*t1;
    i1 = integral(@(t)I1(t).*exp(Ai(2,2)*(t2-t)),t1,t2);
    i2 = integral(@(t)I2(t).*exp(Ai(2,2)*(t2-t)),t1,t2);
    num = (-dun*t2-un+dun*t1);
    deno = (t2^2-2*t1*t2+t1^2);    
    
    Ki = Ai(1,1) + i1/deno;    
    zi = Ai(1,2)*yn(2,end)*exp(Ai(2,2)*(tt(i)-tn(end))) + Ai(1,2)*Ai(2,1)*(num/deno*i1+i2);
    
    M = 1-a_rk(i,i)*dt*Ki;
    rhs = yn(1,end) + a_rk(i,i)*dt*zi;
    for j=1:i-1,
        rhs = rhs + dt*a_rk(i,j)*k(:,j);
    end
    Y(1,i) = M\rhs;
    k(i) = Ki*Y(1,i) + zi;
end
Y(2,end) = yn(2)*exp(Ai(2,2)*(tt(end)-tn)) + Ai(2,1)*(Y(1,end)*i1/deno+num/deno*i1+i2);

out=Y(:,end);