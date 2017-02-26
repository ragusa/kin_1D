function rho_test

clear;clc;close all;

syms x t
syms S_p(x,t) phi(x,t) C(x,t) S_c(x,t)
syms varphi(x,t) amp(t) dvarphi(x,t)

syms f(t) a(x,t) A(t)
f(t) = (t+1);
a(x,t) = x*(1-x)*(x+t);
A(t) = int(a(x,t)*a(x,0),x,0,1);
phi(x,t) = f(t)*a(x,t)/A(t);
varphi(x,t) = a(x,t)/A(t);
amp(t) = f(t);

iv = 1e-3; Sa = 1.0; nfSf = 1.1; cdiff = 1.0; b = 0.0;
S_p(x,t) = iv*diff(phi,t) + (Sa-nfSf*(1-b))*phi - cdiff*diff(diff(phi,x),x);

phi_exact = matlabFunction(phi);
shape_exact = matlabFunction(varphi);
dshape_exact = matlabFunction(diff(varphi,t));
p_exact = matlabFunction(amp);
source_phi = matlabFunction(S_p);
clear x t S_p phi C S_c 
clear varphi(x,t) amp(t) dvarphi(x,t)
clear f(t) a(x,t) A(t)

x = [0 1/3 2/3 1]';
tend = 1.0;

pinit = 1.0;
phi_adj = phi_exact(x,0);
K0 = 0.105;

PL = zeros(4);
PL(2,2) = -10.7614285714286;
PL(3,2) =  7.42017857142857;
PL(2,3) =  7.42017857142857;
PL(3,3) = -10.7614285714286;

rho_exact = @(t) (phi_adj' * PL * shape_exact(x,t)) / K0;
q = @(t) (phi_adj' * assemble_source(source_phi,t))/K0;

n_micro = 10;
n_macro = [10 20 40 80 500 1000];% 5000 1e4];
ntimes=n_macro*n_micro; % 5e4 1e5];% 5e5];% 5e5];

% rho_fun = @(t,told,tnew) rho_exact(t); type = 'Exact \rho(t)';
% rho_fun = @(t,told,tnew) lagrange_fun(rho_exact,told,tnew,1,t); type = 'Linear interpolation of \rho(t) over macro step';
% rho_fun = @(t,told,tnew) lagrange_fun(rho_exact,told,tnew,2,t); type = 'Quadratic lagrange interpolation of \rho(t) over macro step';
% rho_fun = @(t,told,tnew) lagrange_fun(rho_exact,told,tnew,3,t); type = 'Cubic lagrange interpolation of \rho(t) over macro step';
% rho_fun = @(t,told,tnew) (phi_adj' * PL * lagrange_fun(@(tt) shape_exact(x,tt),told,tnew,1,t)) / K0; type = 'Linear interpolation of shape(t) over macro step';
% rho_fun = @(t,told,tnew) (phi_adj' * PL * lagrange_fun(@(tt) shape_exact(x,tt),told,tnew,2,t)) / K0; type = 'Quadratic Lagrange interpolation of shape(t) over macro step';
% rho_fun = @(t,told,tnew) (phi_adj' * PL * lagrange_fun(@(tt) shape_exact(x,tt),told,tnew,3,t)) / K0; type = 'Cubic Lagrange interpolation of shape(t) over macro step';
% rho_fun = @(t,told,tnew) (phi_adj' * PL * hermite_fun(shape_exact(x,told),shape_exact(x,tnew),dshape_exact(x,told),dshape_exact(x,tnew),told,tnew,t)) / K0; type = 'Cubic Hermite interpolation of shape(t) over macro step';
% rho_fun = @(t,told,tnew) (phi_adj' * PL * lagrange_fun(@(tt) shape_exact(x,tt),told-(tnew-told)*1,tnew,2,t)) / K0; type = 'Quadratic lagrange interpolation of shape(t) from previous macro steps';
% rho_fun = @(t,told,tnew) (phi_adj' * PL * lagrange_fun(@(tt) shape_exact(x,tt),told-(tnew-told)*2,tnew,3,t)) / K0; type = 'Cubic lagrange interpolation of shape(t) from previous macro steps';


do_bdf2=true;
do_bdf3=true;


time_integration='IE';
switch time_integration
    case 'IE'
        rk.nstages=1;
        rk.g=1.0;
        rk.a=[rk.g];
        rk.c=sum(rk.a,2);
        rk.b=rk.a(rk.nstages,:);
        rk.order=1;
    case 'sdirk3'
        rk.nstages=3;
        rk.g=0.43586652150845899941601945119356;
        rk.a=[rk.g 0 0; ...
            ((1-rk.g)/2) rk.g 0;...
            (-(6*rk.g^2-16*rk.g+1)/4) ((6*rk.g^2-20*rk.g+5)/4) rk.g];
        rk.c=sum(rk.a,2);
        rk.b=rk.a(rk.nstages,:);
        rk.order=3;
    case 'sdirk3a'
        rk.nstages=3;
        rk.g=0.2113248654051871177454256097490213;
        rk.a=[rk.g, 0, 0 ;...
            rk.g, rk.g, 0;...
            rk.g, 0.5773502691896257645091487805019573, rk.g ];
        rk.c=sum(rk.a,2);
        rk.b=rk.a(rk.nstages,:);
        rk.order=3;
    case 'sdirk2a'
        rk.nstages=2;
        rk.g=0.2928932188134524755991556378951510;
        rk.a=[rk.g, 0 ;...
            1-rk.g, rk.g];
        rk.c=sum(rk.a,2);
        rk.b=rk.a(rk.nstages,:);
        rk.order=2;
    case 'sdirk2b'
        rk.nstages=2;
        rk.g=1.707106781186547524400844362104849;
        rk.a=[rk.g, 0 ;...
            1-rk.g, rk.g];
        rk.c=sum(rk.a,2);
        rk.b=rk.a(rk.nstages,:);
        rk.order=2;
    case 'sdirk4b'
        rk.nstages=5;
        rk.g=0.25;
        rk.a=[rk.g, 0, 0, 0, 0 ;...
            0.5,rk.g, 0, 0, 0;...
            0.34, -0.04, rk.g, 0, 0;...
            0.2727941176470588235294117647058824, -0.05036764705882352941176470588235294, 0.02757352941176470588235294117647059, rk.g, 0;...
            1.041666666666666666666666666666667,-1.020833333333333333333333333333333,7.8125,-7.083333333333333333333333333333333, rk.g ];
        rk.c=sum(rk.a,2);
        rk.b=rk.a(rk.nstages,:);
        rk.order=4;
    case 'sdirk4a'
        rk.nstages=5;
        rk.g=0.2666666666666666666666666666666667;
        rk.a=[rk.g, 0, 0, 0, 0 ;...
            0.5,rk.g, 0, 0, 0;...
            0.3541539528432732316227461858529820, -0.05415395284327323162274618585298197, rk.g, 0, 0;...
            0.08515494131138652076337791881433756, -0.06484332287891555171683963466229754, 0.07915325296404206392428857585141242, rk.g, 0;...
            2.100115700566932777970612055999074, -0.7677800284445976813343102185062276, 2.399816361080026398094746205273880, -2.998818699869028161397714709433394, rk.g ];
        rk.c=sum(rk.a,2);
        rk.b=rk.a(rk.nstages,:);
        rk.order=4;
    otherwise
        error('this sdirk is not available');
end


for iconv=1:length(ntimes)
    dt=tend/ntimes(iconv);
    tn=0;
    pn=pinit;
    fprintf('%d/%d\n',iconv,length(ntimes));
    y=zeros(rk.nstages,1); SS=y; RR=y;
    for it=1:n_macro(iconv)
        rho = @(t) rho_fun(t,tn,tn+dt*n_micro);
        for jt=1:n_micro
            for i=1:rk.nstages
                ti = tn + rk.c(i)*dt;
                RR(i)=rho(ti);
                SS(i)=q(ti);
                rhs = pn + dt*rk.a(i,i)*SS(i);
                for j=1:i-1
                    tj = tn + rk.c(j)*dt;
                    rhs = rhs + dt*rk.a(i,j)*(RR(j)*y(j)+SS(j));
                end
                y(i) = (1 -rk.a(i,i)*dt*RR(i)) \ rhs;
            end
            pn = y(end);
            tn = tn + dt;
        end
    end
    err(iconv)=abs(pn-p_exact(tend))/p_exact(tend);
end

loglog(tend./ntimes, err, '+-' )
p = polyfit(log10(tend./ntimes),log10(err),1);
% slope_err = polyfit(log10(tend./ntimes), log10(err),1)
% reference line y=Cx^s log(y) = log(C) + s log(x)
% we pick one value of (x,y) to get C
% the multiplier in front of C is used to shift the ref. line
% nn=1; %
% nn=length(err);
% s=rk.order; C=0.5*err(nn)/(tend/ntimes(nn))^s;
% y=C*(tend./ntimes).^s;
% hold all;loglog(tend./ntimes, y, 'r-')
leg=[time_integration ', slope = ' num2str(p(1))]; 
% leg=char(leg,sprintf('slope %d',s));
hold on;

if do_bdf2
    for iconv=1:length(ntimes)
        dt=tend/ntimes(iconv);
        tn=0;
        pn=pinit;
        fprintf('%d/%d\n',iconv,length(ntimes));
        polder=pn;
        % first time step: CN
        %         for it=1:ntimes(iconv)
        rho = @(t) rho_fun(t,tn,tn+dt*n_micro);
        pn = (1-dt/2*rho(tn+dt))\ ( (1+dt/2*rho(tn))*pn + dt/2*(q(tn)+q(tn+dt)));
        %         end
        % other time steps: BDF2
        for it=1:n_macro(iconv)
            rho = @(t) rho_fun(t,tn,tn+dt*n_micro);
            for jt=1:n_micro
                if it>1 || jt>1
                    pnew = ( 1 - 2/3*dt*rho(tn+dt)) \ ( 1/3*(4*pn-polder) +2/3*dt*q(tn+dt));
                    polder=pn;
                    pn=pnew;
                end
                tn = tn + dt;
            end
        end
        err(iconv)=abs(pn-p_exact(tend))/p_exact(tend);
    end
    
    loglog(tend./ntimes, err, 'o-' )
    p = polyfit(log10(tend./ntimes), log10(err),1);
    % reference line y=Cx^s log(y) = log(C) + s log(x)
    % we pick one value of (x,y) to get C
    % the multiplier in front of C is used to shift the ref. line
%     nn=1; %
%     nn=length(err);
%     s=2; C=0.5*err(nn)/(tend/ntimes(nn))^s;
%     y=C*(tend./ntimes).^s;
%     hold all;loglog(tend./ntimes, y, 'r-')
    leg=char(leg,['BDF2, slope = ' num2str(p(1))]); 
%     leg=char(leg,sprintf('slope %d',2));
end

if do_bdf3
    for iconv=1:length(ntimes)
        dt=tend/ntimes(iconv);
        tn=0;
        pn=pinit;
        fprintf('%d/%d\n',iconv,length(ntimes));
        % first time step: CN
        rho = @(t) rho_fun(t,tn,tn+dt*n_micro);
        pn = (1-dt/2*rho(tn+dt))\ ( (1+dt/2*rho(tn))*pn + dt/2*(q(tn)+q(tn+dt)));
        % second time step: BDF2
        polder=pinit;
        pnew = ( 1 -2/3*dt*rho(tn+2*dt)) \ ( 1/3*(4*pn-polder) +2/3*dt*q(tn+2*dt));
        poldest=pinit;
        polder=pn;
        pn=pnew;
        % other time steps: BDF3
        for it=1:n_macro(iconv)
            rho = @(t) rho_fun(t,tn,tn+dt*n_micro);
            for jt=1:n_micro
                if it>1 || jt>2
                    pnew = ( 1 -6/11*dt*rho(tn+dt)) \ ( 1/11*(18*pn-9*polder+2*poldest) +6/11*dt*q(tn+dt));
                    poldest=polder;
                    polder=pn;
                    pn=pnew;
                end
                tn = tn + dt;
            end
        end
        err(iconv)=abs(pn-p_exact(tend))/p_exact(tend);
    end
    
    loglog(tend./ntimes, err, 's-' )
    p = polyfit(log10(tend./ntimes), log10(err),1);
    % reference line y=Cx^s log(y) = log(C) + s log(x)
    % we pick one value of (x,y) to get C
    % the multiplier in front of C is used to shift the ref. line
%     nn=1; %
%     nn=length(err);
%     s=3; C=0.5*err(nn)/(tend/ntimes(nn))^s;
%     y=C*(tend./ntimes).^s;
%     hold all;loglog(tend./ntimes, y, 'r-')
    leg=char(leg,['BDF3, slope = ' num2str(p(1))]); 
%     leg=char(leg,sprintf('slope %d',3));
    
end
   
legend(leg,'Location','Best')
grid on;
xlabel('\Delta t_m_i_c_r_o')
ylabel('Error')
title(type)

return
end

function val = lagrange_fun(fun,told,tnew,order,t)
dt = (tnew-told)/order;
tt = told:dt:tnew;
val = 0.0;
for i=1:order+1
    P = fun(tt(i));
    for j=1:order+1
        if i~=j
            P = P.*(t-tt(j))/(tt(i)-tt(j));
        end
    end
    val = val + P;
end
return
end

function val = hermite_fun(fun_beg,fun_end,dfun_beg,dfun_end,t1,t2,t)
mat=[ t1^3 t1^2 t1 1; t2^3 t2^2 t2 1;  3*t1^2 2*t1 1 0; 3*t2^2 2*t2 1 0];
val = zeros(length(fun_beg),1);
for i=1:length(fun_beg)
    rhs = [ fun_beg(i) ;fun_end(i) ; dfun_beg(i); dfun_end(i)];
    w = mat\rhs;
    val(i) = [t^3 t^2 t 1]*w;
end
return
end