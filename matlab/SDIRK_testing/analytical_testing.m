clc; clear; close all; hold all

t_end = 1.0;
min_N = 6;
max_N = 12;
N = 2.^(min_N:max_N);
error_bf = zeros(1,length(N)); error_an = error_bf;

gamma = 0.43586652150846;
a = 1/2*(1-gamma);
b = .25*(-6*gamma^2+16*gamma-1);
c = .25*(6*gamma^2-20*gamma+5);
a_rk = [gamma 0 0; a gamma 0; b c gamma];
b_rk = [b c gamma];
c_rk = [gamma ; a+gamma ; 1];

A = [4 1; 2 -1];
A_fun = @(t) [4 1; 2 -1];
X0 = [1; 1];

options = odeset('RelTol',1e-12,'AbsTol',1e-12,'InitialStep',1e-8,'Stats','on','MaxStep',1e-1,'Stats','off');
[t_ref,X_ref] = ode15s(@Ax_function,[0 t_end], X0,options);
figure(1)
plot(t_ref,X_ref(:,1));
leg = char('Reference');
legend(leg,'Location','Best');


for i=1:length(N)
    dt = t_end/N(i);
    t = 0:dt:t_end;
    X_bf = zeros(2,length(t)); X_bf(:,1) = X0;
    X_an = zeros(2,length(t)); X_an(:,1) = X0;
    
    for n=1:N(i)
        X_bf(:,n+1) = SDIRK(X_bf(:,n),t(n),A_fun,a_rk,c_rk,dt);
        X_an(:,n+1) = SDIRK_an(X_an(:,n),t(n),A_fun,a_rk,c_rk,dt);
    end
    
    figure(1)
    plot(t,X_bf(1,:),'color',rand(1,3))
    plot(t,X_an(1,:),'color',rand(1,3))
    leg = char(leg,strcat('Brute Force N = ',num2str(N(i))),strcat('Analytical N = ',num2str(N(i))));
    legend(leg,'Location','Best');
    
    error_bf(i) = abs(X_bf(1,end)-X_ref(end,1))/X_ref(end,1);
    error_an(i) = abs(X_an(1,end)-X_ref(end,1))/X_ref(end,1);
end

p_bf = polyfit(log10(t_end./N),log10(error_bf),1);
p_an = polyfit(log10(t_end./N),log10(error_an),1);

figure(2)
loglog(t_end./N,error_bf,'o-',t_end./N,error_an,'o-')
legend(strcat('Brute Force, Slope = ',num2str(p_bf(1))),strcat('Analytical, Slope = ',num2str(p_an(1))));
hold off
