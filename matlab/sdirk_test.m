clear all; close all; clc;

bf=@(x)4*x.*(1-x);
a=integral(@(x)bf(x).*bf(x),0,1)
b=integral(@(x)bf(x),0,1)
g=integral(@(x)(4*(1-2*x)).^2,0,1)
D=1;
v=1;
nsf=1.1;sa=1;
p=2;

A=-D*g+(nsf-sa)*a;
IV=a/v;
S =@(t) a/4/v*p*(1+t).^(p-1) + (1+t).^p*((sa-nsf)*a/4+2*D*b);
ex=@(t)(1+t).^p/4;

rk.nstages=3;
rk.g=0.43586652150845899941601945119356;
rk.a=[rk.g 0 0; ...
    ((1-rk.g)/2) rk.g 0;...
    (-(6*rk.g^2-16*rk.g+1)/4) ((6*rk.g^2-20*rk.g+5)/4) rk.g];
rk.c=sum(rk.a,2);
rk.b=rk.a(rk.nstages,:);

yold=0.25;

tend=.1;

ntimes=[10 20 40];
for iconv=1:length(ntimes)
    dt=tend/ntimes(iconv); 
    tn=0;
    yn=yold;
    for it=1:ntimes(iconv)
        for i=1:rk.nstages
            ti = tn + rk.c(i)*dt;
            rhs = yn + dt*rk.a(i,i)*S(ti);
            for j=1:i-1
                rhs = rhs + dt*rk.a(i,j)*(A*y(j)+S(tn+dt*rk.c(j)));
            end
            y(i) = (IV -rk.a(i,i)*dt*A) \ rhs;
        end
        yn = y(3);
        tn = tn + dt;
    end
    err(iconv)=abs(yn-ex(tend))
end

plot( log10(tend./ntimes), log10(err) )