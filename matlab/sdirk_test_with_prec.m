clear all; close all; clc;
% exact solution for flux: (1-x)*x*(1+t)^p = bf2(x)/4*(1+t)^p
% exact solution for precusors: (1-x)*x*(1+t)^p = bf2(x)/4*(1+t)^p
p=4;
% pb data
D=1; v=1e4;
beta=0.001; nsf=1.1;sa=1;
lambda=0.1;

% first basis function
bf1=@(x)(1-x).*(0.5-x)*2;
% last basis function
bf3=@(x)x.*(x-0.5)*2;
% middle basis function
bf2=@(x)4*x.*(1-x);

% fem entry for second basis function
a=integral(@(x)bf2(x).*bf2(x),0,1);
b=integral(@(x)bf2(x),0,1);
g=integral(@(x)(4*(1-2*x)).^2,0,1);
A=-D*g+(nsf*(1-beta)-sa)*a;
IV=a/v;
% source term for the flux equation after fem stuff
S_flx =@(t) a/4/v*p*(1+t).^(p-1) + (1+t).^p*((sa-nsf*(1-beta))*a/4 + 2*D*b) + lambda*a/4*(1+t).^p;
% exact solution for the flux at x=0.5 ==> bf2(0.5)/4*(1+t)^p = (1+t)^p /4
ex_flux=@(t)(1+t).^p/4;

% source term for the precursors after fem stuff (precursors that are not coupled to flux, 1 and 3)
b1b2=integral(@(x)bf1(x).*bf2(x),0,1);
S_p1 = @(t) p*(1+t).^(p-1)*b1b2/4 + (1+t).^p*lambda*b1b2/4;
S_p2 = @(t) p*(1+t).^(p-1)*a/4    + (1+t).^p*lambda*a/4    -beta*nsf*a/4*(1+t).^p;
% no need for 3 equation because b1b2=b3b2
b3b2=integral(@(x)bf3(x).*bf2(x),0,1);
ex_p1=@(t)(1+t).^p *b1b2/4;
ex_p2=@(t)(1+t).^p *a/4;

yinit=[bf2(0.5)/4  b1b2/4 a/4]';
tend=5;

TD_matrix = [ IV 0 0 ; 0 1 0 ; 0 0 1];
big_A = [ A lambda lambda ; beta*nsf*b1b2/4 -lambda 0; beta*nsf*a/4 0 -lambda]; 
S_vec = @(t) [ S_flx(t) ; S_p1(t) ; S_p2(t) ];

do_bdf2=true;
do_bdf3=true;

ntimes=[10 20 40 80 500 1000 5000 1e4]; % 5e4 1e5];% 5e5];% 5e5];

time_integration='sdirk4a';
switch time_integration
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
    yn=yinit;
    fprintf('%d/%d\n',iconv,length(ntimes));
    y=zeros(rk.nstages,1); SS=y;
    for it=1:ntimes(iconv)
        for i=1:rk.nstages
            ti = tn + rk.c(i)*dt;
            SS(i)=S(ti);
            rhs = IV*yn + dt*rk.a(i,i)*SS(i);
            for j=1:i-1
                rhs = rhs + dt*rk.a(i,j)*(A*y(j)+SS(j));
            end
            y(i) = (IV -rk.a(i,i)*dt*A) \ rhs;
        end
        yn = y(end);
        tn = tn + dt;
    end
    err(iconv)=abs(yn-ex(tend));
end

plot( log10(tend./ntimes), log10(err), '+-' )
% slope_err = polyfit(log10(tend./ntimes), log10(err),1)
% reference line y=Cx^s log(y) = log(C) + s log(x)
% we pick one value of (x,y) to get C
% the multiplier in front of C is used to shift the ref. line
nn=1; %
nn=length(err);
s=rk.order; C=0.5*err(nn)/(tend/ntimes(nn))^s;
y=C*(tend./ntimes).^s;
hold all;plot(log10(tend./ntimes), log10(y), 'r-')
leg=char(time_integration); 
leg=char(leg,sprintf('slope %d',s)); 


if do_bdf2
    leg=char(leg,'BDF2'); 
    leg=char(leg,sprintf('slope %d',2));
    for iconv=1:length(ntimes)
        dt=tend/ntimes(iconv);
        tn=0;
        yn=yinit;
        fprintf('%d/%d\n',iconv,length(ntimes));
        y=zeros(rk.nstages,1); SS=y;
        yolder=yn;
        % first time step: CN
        %         for it=1:ntimes(iconv)
        yn = (IV-dt/2*A)\ ( (IV+dt/2*A)*yn + dt/2*(S(tn)+S(tn+dt)));
        tn = tn + dt;
        %         end
        % other time steps: BDF2
        for it=2:ntimes(iconv)
            ynew = ( IV -2/3*dt*A) \ ( IV/3*(4*yn-yolder) +2/3*dt*S(tn+dt));
            tn = tn + dt;
            yolder=yn;
            yn=ynew;
        end
        err(iconv)=abs(yn-ex(tend));
    end
    
    plot( log10(tend./ntimes), log10(err), 'o-' )
    % slope_err = polyfit(log10(tend./ntimes), log10(err),1)
    % reference line y=Cx^s log(y) = log(C) + s log(x)
    % we pick one value of (x,y) to get C
    % the multiplier in front of C is used to shift the ref. line
    nn=1; %
    nn=length(err);
    s=2; C=0.5*err(nn)/(tend/ntimes(nn))^s;
    y=C*(tend./ntimes).^s;
    hold all;plot(log10(tend./ntimes), log10(y), 'r-')
    
end



if do_bdf3
    leg=char(leg,'BDF3'); 
    leg=char(leg,sprintf('slope %d',3));
    for iconv=1:length(ntimes)
        dt=tend/ntimes(iconv);
        tn=0;
        yn=yinit;
        fprintf('%d/%d\n',iconv,length(ntimes));
        y=zeros(rk.nstages,1); SS=y;
        % first time step: CN
        yn = (IV-dt/2*A)\ ( (IV+dt/2*A)*yn + dt/2*(S(tn)+S(tn+dt)));
        tn = tn + dt;
        % second time step: BDF2
        yolder=yinit;
        ynew = ( IV -2/3*dt*A) \ ( IV/3*(4*yn-yolder) +2/3*dt*S(tn+dt));
        tn = tn + dt;
        yoldest=yinit;
        yolder=yn;
        yn=ynew;
        % other time steps: BDF3
        for it=3:ntimes(iconv)
            ynew = ( IV -6/11*dt*A) \ ( IV/11*(18*yn-9*yolder+2*yoldest) +6/11*dt*S(tn+dt));
            tn = tn + dt;
            yoldest=yolder;
            yolder=yn;
            yn=ynew;
        end
        err(iconv)=abs(yn-ex(tend));
    end
    
    plot( log10(tend./ntimes), log10(err), 's-' )
    % slope_err = polyfit(log10(tend./ntimes), log10(err),1)
    % reference line y=Cx^s log(y) = log(C) + s log(x)
    % we pick one value of (x,y) to get C
    % the multiplier in front of C is used to shift the ref. line
    nn=1; %
    nn=length(err);
    s=3; C=0.5*err(nn)/(tend/ntimes(nn))^s;
    y=C*(tend./ntimes).^s;
    hold all;plot(log10(tend./ntimes), log10(y), 'r-')
    
end

legend(leg)