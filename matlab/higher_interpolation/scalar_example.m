clear all; close all; clc;

% we solve du/dt = a(t).u(t) + q(t) with u(t=0)=u0
% a(t) is an arbitrary function of time
% when employing a user-chosen exact solution, u=ue, then this
% automatically sets q(t)= du/dt-a.u. otherwise, q(t)=0

% end of simulation time
tend=100;

% coefficients for time dependence of exact solution, if selected
coef = [1 1 1 1 1 1]; % coef for t^0, t^1, t^2, t^3, t^4, t^5
coef = [0 0 0 0 1 1]; % coef for t^0, t^1, t^2, t^3, t^4, t^5

% a(t) choice 0(constant) 1(linear) 2(decay exp) 3(incr exp)
a_type = 0;
a_const=0; % used only if a_type=0

% solution type: 0(no mms) 1(linear) 2(quadratic) 3(cubic) 4(qurtic)
% 5(quintic)
solution_type = 4;

% time_discretization 'Crank-Nicholson' or 'SDIRK33' or 'SDIRK54'
time_discretization = 'Crank-Nicholson';
time_discretization = 'SDIRK54';

% sdirk33 constants:
g=0.43586652150845899941601945119356;
A=[g 0 0; ...
    ((1-g)/2) g 0;...
    (-(6*g^2-16*g+1)/4) ((6*g^2-20*g+5)/4) g];
c=sum(A'); b=A(end,:);

% sdirk54 constants:
A=[ 1./4., 0., 0., 0., 0.;...
    1./2., 1./4., 0., 0., 0.;...
    17./50., -1./25., 1./4., 0., 0.;...
    371./1360., -137./2720., 15./544., 1./4., 0.;...
    25./24., -49./48., 125./16., -85./12., 1./4.];
c=sum(A'); b=A(end,:);

% tolerances for odesolvers
rtol = 1e-13; abso = 1e-13;
atol  = abso*ones(length([1]),1);
options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);

% pick function for a(t)
switch a_type
    case 0 % constant
        a = @(t)(0.*t+a_const);
        dadt = @(t)(0.*t);
    case 1 % linear
        a = @(t)(t);
        dadt = @(t)(0.*t+1.);
    case 2 % decaying exponential
        a = @(t)(exp(-t));
        dadt = @(t)(-exp(-t));
    case 3 % increasing exponential
        a = @(t)(exp(t));
        dadt = @(t)(exp(-t));
    otherwise
        error('unknown type for function a(t)');
end

% pick function for q(t)
if solution_type == 0
    % no manufactured solution
    % set q(t)
    q = @(t) (0.*t);
    % impose initial value
    u0 = 1;
    % solve the ODE to get the final time reference value
    [tt,uu]=ode15s(@scalar_ssres,[0 tend],u0,options,a,q);
    % extract reference value
    u_ref = uu(end);
    % plot
    plot(tt,uu); hold all;
else
    switch solution_type
        
        case 1 % u_exact is linear
            exact = @(t)(coef(1)+coef(2)*t);
            dexactdt = @(t)(coef(2));
            
        case 2 % u_exact is quadratic
            exact = @(t)(coef(1)+coef(2)*t+coef(3)*t.^2);
            dexactdt = @(t)(coef(2)+2*coef(3)*t);
            
        case 3 % u_exact is cubic
            exact = @(t)(coef(1)+coef(2)*t+coef(3)*t.^2+coef(4)*t.^3);
            dexactdt = @(t)(coef(2)+2*coef(3)*t+3*coef(4)*t.^2);
            
        case 4 % u_exact is quartic
            exact = @(t)(coef(1)+coef(2)*t+coef(3)*t.^2+coef(4)*t.^3+coef(5)*t.^4);
            dexactdt = @(t)(coef(2)+2*coef(3)*t+3*coef(4)*t.^2+4*coef(5)*t.^3);
            
        case 5 % u_exact is quintic
            exact = @(t)(coef(1)+coef(2)*t+coef(3)*t.^2+coef(4)*t.^3+coef(5)*t.^4+coef(6)*t.^5);
            dexactdt = @(t)(coef(2)+2*coef(3)*t+3*coef(4)*t.^2+4*coef(5)*t.^3+5*coef(6)*t.^4);
            
        otherwise
            error('unknown type for solution');
    end
    % set q(t)
    q = @(t) (dexactdt(t) - a(t)*exact(t));
    % reference value
    u_ref = exact(tend);
    % plot
    tt=linspace(0,tend);
    plot(tt,exact(tt)); hold all;
    % impose initial value
    u0 = exact(0);
    % solve the ODE to get the final time reference value
    [tt,uu]=ode15s(@scalar_ssres,[0 tend],u0,options,a,q);
    % plot
    plot(tt,uu);
    fprintf('u_ref = %15.12e, u_ode = %15.12e, errpr = %g \n',u_ref,uu(end),u_ref-uu(end));
end
% u_ref = uu(end);

%%% convergence study
n_runs=12;
n_steps = 2.^linspace(0,n_runs-1,n_runs);

% n_steps=1;

for i_run=1:n_runs
    
    nbr_steps = n_steps(i_run);
    dt = tend/nbr_steps;
    
    sol=zeros(nbr_steps+1,1);
    sol(1)=u0;
    
    for it=1:nbr_steps
        if strcmpi(time_discretization,'Crank-Nicholson')
            time0 = (it-1)*dt;
            time1 = time0 + dt;
            % CN: (unew - uold)/dt = 0.5 ( a(old)u(old)+q(old) + a(new)u(new)+q(new) )
            deno = 1-dt/2*a(time1);
            sol(it+1) = ( sol(it) + dt/2*(scalar_ssres(time0,sol(it),a,q) + q(time1)) )/deno;
        elseif strcmpi(time_discretization(1:5),'SDIRK')
            % number of stages
            n_stages = str2double(time_discretization(6));
            Y=zeros(n_stages,1); ts=Y;
            % time at the beginning of the time step
            time0 = (it-1)*dt;
            
            % Yi = yn + dt sum_j { A_ij f(tj, Yj) }
            for i=1:n_stages
                % stage time
                ts(i) = time0 + c(i)*dt;
                % denominator
                deno = 1 - dt*A(i,i)*a(ts(i));
                % numerator
                nume = sol(it) + dt*A(i,i)*q(ts(i)) ;
                for j=1:i-1
                    nume = nume + dt*A(i,j)*scalar_ssres(ts(j),Y(j),a,q);
                end
                Y(i) = nume/deno;
                %                 % stage 1
                %                 Y1 = ( sol(it) + dt*A(1,1)*q(t1) )/deno;
                %                 % stage 2
                %                 t2 = time0 + c(2)*dt;
                %                 deno = 1 - dt*A(2,2)*a(t2);
                %                 Y2 = ( sol(it) + dt*A(2,1)*scalar_ssres(t1,Y1,a,q) + dt*A(2,2)*q(t2) )/deno;
                %                 % stage 3
                %                 t3 = time0 + c(3)*dt;
                %                 deno = 1 - dt*A(3,3)*a(t3);
                %                 Y3 = ( sol(it) + dt*A(3,1)*scalar_ssres(t1,Y1,a,q) + dt*A(3,2)*scalar_ssres(t2,Y2,a,q) + dt*A(3,3)*q(t3) )/deno;
            end
            % update
            sol(it+1) = Y(n_stages);
        else
            error('unknown time discretization');
        end
    end
    plot(linspace(0,tend,nbr_steps+1),sol);
    conv_sol(i_run) = abs(sol(end)-u_ref)/u_ref;
    conv_dt(i_run) = dt;
    
end

figure(2)
loglog(conv_dt,conv_sol,'+-'); hold all;
% error = C dt^p. Compute C for the smallest dt value
if strcmpi(time_discretization,'Crank-Nicholson')
    porder=2;
elseif strcmpi(time_discretization,'SDIRK33')
    porder=3;
elseif strcmpi(time_discretization,'SDIRK54')
    porder=4;
end
C = conv_sol(end)/conv_dt(end)^porder;
loglog(conv_dt,C*conv_dt.^porder);


[conv_dt' conv_sol']