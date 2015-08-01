function kinetics_driver2

clear all;
close all; clc;

global dat npar

% verbose/output parameters
testing = false;
console_print = true;
plot_transient_figure = true;
plot_power_figure = false;
make_movie = false;

% one of the two choices for applying BC
npar.set_bc_last=true;

% select problem
pbID=2; refinements=20;
problem_init(pbID,refinements);
npar.mf_sol = true;

% compute fundamental eigenmode
curr_time=0;
[phi0,keff]=steady_state_eigenproblem(curr_time);
% initialize kinetic values
C0 = kinetics_init(phi0,curr_time);
if npar.mf_sol
    phi0 = (npar.x_dofs/dat.width.*(1-npar.x_dofs/dat.width))';
    C0 = (npar.x_dofs/dat.width.*(1-npar.x_dofs/dat.width))';
    npar.phi_adj = phi0;
end

% initial solution vector
u=[phi0;C0];
% save a copy of it
u0=u;

phi_adjoint = npar.phi_adj;
IV   = assemble_mass(     dat.inv_vel ,curr_time);
npar.IVel = IV;

if plot_transient_figure
    plot(npar.x_dofs,phi0);
end
if console_print
    fprintf('%10.8g \n',keff);
end

% time steping data
dt=0.001;
ntimes=2000; % 150*2;
iqs_factor=1;


if make_movie
    %# figure
    figure, set(gcf, 'Color','white')
    axis([0 400 0 0.65]);
    set(gca, 'nextplot','replacechildren', 'Visible','off');
    %# preallocate
    mov(1:ntimes) = struct('cdata',[], 'colormap',[]);
end

% save flux for tests
phi_save=zeros(length(phi0),ntimes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% brute force discretization of
%%%   the TD neutron diffusion eq and precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dobr = true;
if dobr
amplitude_norm=1;

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    
    % solve time-dependent diffusion for flux
    u = solve_TD_diffusion(u,dt,time_end);
    
    % compute end time power for plotting
    dat.Ptot(it+1) = compute_power(dat.nusigf,time_end,u(1:npar.n));
    
    % ratio of <u*,IVu> to its initial value
    amplitude_norm(it+1) = (phi_adjoint'*IV*u(1:npar.n))/npar.K0;    
    
    % plot/movie
    if plot_transient_figure
        figure(1);
        if npar.mf_sol
            plot(npar.x_dofs,u(1:npar.n),'b-',npar.x_dofs,npar.x_dofs/dat.width.*(1-npar.x_dofs/dat.width)*(time_end+1),'ko');drawnow;
        else
            plot(npar.x_dofs,u(1:npar.n));drawnow;
        end
        if make_movie, mov(it) = getframe(gca); end
    end

    
    % save flux for usage in PRKE exact for testing purposes
    phi_save(:,it)=u(1:npar.n); %/Pnorm(it+1);
    
end
% make movie
if plot_transient_figure && make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, 'PbID10_v2.avi', 'compression','None', 'fps',1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS IQS IQS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doiqs = false;
if doiqs
fig=2;
theta_log=false;
[amplitude_norm_iqs1] = IQS_driver(phi0,C0,dt,ntimes,iqs_factor, plot_transient_figure,theta_log,fig);
end
% fig=3;
% theta_log=true;
% [amplitude_norm_iqs2] = IQS_driver(phi0,C0,dt,ntimes,iqs_factor, plot_transient_figure,theta_log,fig);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% standard PRKE with initial shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doprke = false;
if doprke
shape0=phi0;
[~,beff_MGT]=compute_prke_parameters(0.,shape0);
X=[1;beff_MGT/dat.lambda];
Pnorm_prke(1)=X(1);

% loop over time steps
for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    % solve prke
    X =  solve_prke(X,dt,time_end,shape0);
    % store power level for plotting
    Pnorm_prke(it+1)=X(1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prke using exact shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% shape0=phi0;
% [~,beff_MGT]=compute_prke_parameters(0.,shape0);
% X=[1;beff_MGT/dat.lambda];
% Pnorm_prkeEX(1)=X(1);
% IV = assemble_mass(dat.inv_vel ,0.);
% 
% % loop over time steps
% for it=1:ntimes
%     time_end=it*dt;
%     if console_print, fprintf('time end = %g \n',time_end); end
%     % compute shape from saved flux
%     shape_curr = phi_save(:,it) / (phi_adjoint'*IV*phi_save(:,it)) * npar.K0;
%     % solve prke
%     X =  solve_prke(X,dt,time_end,shape_curr);
%     % store power level for plotting
%     Pnorm_prkeEX(it+1)=X(1);
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% prke using exact QS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% shape0=phi0;
% [~,beff_MGT]=compute_prke_parameters(0.,shape0);
% X=[1;beff_MGT/dat.lambda];
% Pnorm_prkeQS(1)=X(1);
% IV = assemble_mass(dat.inv_vel ,0.);
% 
% % loop over time steps
% for it=1:ntimes
%     time_end=it*dt;
%     if console_print, fprintf('time end = %g \n',time_end); end
%     % compute shape for steady state
%     [shape_curr,~]=steady_state_eigenproblem(time_end);
%     % solve prke
%     X =  solve_prke(X,dt,time_end,shape_curr);
%     % store power level for plotting
%     Pnorm_prkeQS(it+1)=X(1);
% end

if plot_power_figure
    figure(4); hold all;
    t_=linspace(0,dt*ntimes,ntimes+1);
    ti=linspace(0,dt*ntimes,ntimes/iqs_factor+1);
      plot(t_,amplitude_norm,'+-');  leg=char('space-time');
%     plot(t_,Pnorm_prkeEX,'x-');    leg=char(leg,'PRKE exact');
     plot(t_,Pnorm_prke,'ro-');     leg=char(leg,'PRKE');
    plot(ti,amplitude_norm_iqs1,'gx-');  leg=char(leg,'PRKE IQS');
%     plot(ti,amplitude_norm_iqs2,'kx-');  leg=char(leg,'PRKE IQS w/ \theta');
    legend(leg,'Location','Best')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end