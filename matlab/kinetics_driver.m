function kinetics_driver

clear all;
close all; clc;

global dat npar

% verbose/output parameters
testing = true;
console_print = true;
plot_transient_figure = true;
plot_power_figure = true;
make_movie = false;

% one of the two choices for applying BC
npar.set_bc_last=true;

% select problem
pbID=10; refinements=1;
problem_init(pbID,refinements);

% compute fundamental eigenmode
curr_time=0;
[phi0,keff]=steady_state_eigenproblem(curr_time);
if plot_transient_figure
    plot(npar.x_dofs,phi0);
end
if console_print
    fprintf('Initial SS eigenvalue = %10.8g \n',keff);
end

% initialize kinetic values
C0 = kinetics_init(phi0,curr_time);

% initial solution vector
u=[phi0;C0];
% save a copy of it
u0=u;

phi_adjoint = npar.phi_adj;
IV   = assemble_mass(     dat.inv_vel ,curr_time);
npar.IV=IV;

% time steping data
dt=0.005;
ntimes=10/dt; % 150*2;
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

% testing with another weighting function
% npar.phi_adj = ones(length(npar.phi_adj),1);
% npar.phi_adj(1)=0;
% npar.phi_adj(end)=0;

i=0;
% i=i+1; list_runs{i}= 'brute_force';
% i=i+1; list_runs{i}= 'brute_force_an_prec';
i=i+1; list_runs{i}= 'iqs_an_prec';
i=i+1; list_runs{i}= 'iqs';
% i=i+1; list_runs{i}= 'prke_initial_shape';
% i=i+1; list_runs{i}= 'prke_exact_shape';
% i=i+1; list_runs{i}= 'prke_qs_shape';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% brute force discretization of
%%%   the TD neutron diffusion eq and precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'brute_force')
    
amplitude_norm=1;
u=u0;

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    
    % solve time-dependent diffusion for flux
    u = solve_TD_diffusion(u,dt,time_end);
    
    % plot/movie
    if plot_transient_figure
        figure(1);
        plot(npar.x_dofs,u(1:npar.n));drawnow;
        if make_movie, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    dat.Ptot(it+1) = compute_power(dat.nusigf,time_end,u(1:npar.n));
    
    % ratio of <u*,IVu> to its initial value
    amplitude_norm(it+1) = (phi_adjoint'*IV*u(1:npar.n))/npar.K0;
    
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
%%% brute force discretization of
%%%   the TD neutron diffusion eq and ANALYTICAL 
%%%   solution for the precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'brute_force_an_prec')
 
amplitude_norm_an_prec=1;
u=u0;

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    
    % solve time-dependent diffusion for flux
    u = solve_TD_diffusion_an_prec(u,dt,time_end);
    
    % plot/movie
    if plot_transient_figure
        figure(1);
        plot(npar.x_dofs,u(1:npar.n));drawnow;
        if make_movie, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    dat.Ptot_an_prec(it+1) = compute_power(dat.nusigf,time_end,u(1:npar.n));
    
    % ratio of <u*,IVu> to its initial value
    amplitude_norm_an_prec(it+1) = (phi_adjoint'*IV*u(1:npar.n))/npar.K0;
    
    % save flux for usage in PRKE exact for testing purposes
    phi_save_an_prec(:,it)=u(1:npar.n); %/Pnorm(it+1);
    
end
% make movie
if plot_transient_figure && make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, 'PbID10_v2_an_prec.avi', 'compression','None', 'fps',1);
end

% if plot_power_figure
%     figure(3); hold all;
%     t_=linspace(0,dt*ntimes,ntimes+1);
%     plot(t_,amplitude_norm,'+-');  leg=char('space-time');
%     plot(t_,amplitude_norm_an_prec,'+-');  leg=char(leg,'space-time-ANALY');
%     legend(leg,'Location','Best')
% end
% 
% error('done')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS IQS IQS with ANALYTICAL precursors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqs_an_prec')
 
amplitude_norm_iqs2=1;

% initial solution vector
u_shape=[phi0;C0];
[~,beff_MGT]=compute_prke_parameters(0.,phi0);
X=[1;beff_MGT/dat.lambda];
Pnorm_prkeIQS2(1)=X(1);

dt=dt*iqs_factor; ntimes=ntimes/iqs_factor;

n_micro=10;
freq_react=1;

time_prke_iqs=[];
power_prke_iqs=[];

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    
    % solve time-dependent diffusion for flux
    [u_shape,X,t,y] = solve_IQS_diffusion_an_prec(u_shape,X,dt,time_end);
    
    % plot/movie
    if plot_transient_figure
        figure(2)
        plot(npar.x_dofs,u_shape(1:npar.n));drawnow;
        if make_movie, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    dat.Ptot_iqs2(it+1) = compute_power(dat.nusigf,time_end,X(1)*u_shape(1:npar.n));
    
    % ratio of <u*,IVu> to its initial value
    amplitude_norm_iqs2(it+1) = X(1)* (phi_adjoint'*IV*u_shape(1:npar.n))/npar.K0;
    
    % save fine-scale data
    time_prke_iqs=[time_prke_iqs ; t];
    power_prke_iqs=[power_prke_iqs; y];

    
end

% make movie
if plot_transient_figure && make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, 'PbID10_v2_iqs.avi', 'compression','None', 'fps',1);
end
% undo time step increase
dt=dt/iqs_factor; ntimes=ntimes*iqs_factor;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS IQS IQS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqs')
 
amplitude_norm_iqs=1;

% initial solution vector
u_shape=[phi0;C0];
[~,beff_MGT]=compute_prke_parameters(0.,phi0);
X=[1;beff_MGT/dat.lambda];
Pnorm_prkeIQS(1)=X(1);

dt=dt*iqs_factor; ntimes=ntimes/iqs_factor;

n_micro=10;
freq_react=1;

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    
    % solve time-dependent diffusion for flux
    [u_shape,X] = solve_IQS_diffusion(u_shape,X,dt,time_end,n_micro,freq_react);
    
    % plot/movie
    if plot_transient_figure
        figure(2)
        plot(npar.x_dofs,u_shape(1:npar.n));drawnow;
        if make_movie, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    dat.Ptot_iqs(it+1) = compute_power(dat.nusigf,time_end,X(1)*u_shape(1:npar.n));
    
    % ratio of <u*,IVu> to its initial value
    amplitude_norm_iqs(it+1) = X(1)* (phi_adjoint'*IV*u_shape(1:npar.n))/npar.K0;
    
end

% make movie
if plot_transient_figure && make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, 'PbID10_v2_iqs.avi', 'compression','None', 'fps',1);
end
% undo time step increase
dt=dt/iqs_factor; ntimes=ntimes*iqs_factor;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% standard PRKE with initial shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'prke_initial_shape') 

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
if should_I_run_this(list_runs,'prke_exact_shape')
 
shape0=phi0;
[~,beff_MGT]=compute_prke_parameters(0.,shape0);
X=[1;beff_MGT/dat.lambda];
Pnorm_prkeEX(1)=X(1);
IV = assemble_mass(dat.inv_vel ,0.);

% loop over time steps
for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    % compute shape from saved flux
    shape_curr = phi_save(:,it) / (phi_adjoint'*IV*phi_save(:,it)) * npar.K0;
    % solve prke
    X =  solve_prke(X,dt,time_end,shape_curr);
    % store power level for plotting
    Pnorm_prkeEX(it+1)=X(1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prke using exact QS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'prke_qs_shape')

shape0=phi0;
[~,beff_MGT]=compute_prke_parameters(0.,shape0);
X=[1;beff_MGT/dat.lambda];
Pnorm_prkeQS(1)=X(1);
IV = assemble_mass(dat.inv_vel ,0.);

% loop over time steps
for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    % compute shape for steady state
    [shape_curr,~]=steady_state_eigenproblem(time_end);
    % solve prke
    X =  solve_prke(X,dt,time_end,shape_curr);
    % store power level for plotting
    Pnorm_prkeQS(it+1)=X(1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_power_figure
    figure(3); hold all;
    t_=linspace(0,dt*ntimes,ntimes+1);
    ti=linspace(0,dt*ntimes,ntimes/iqs_factor+1);
    has_leg=false;
    
    if should_I_run_this(list_runs,'brute_force')
        plot(t_,amplitude_norm,'+-');
        if ~has_leg
            leg=char('space-time');
            has_leg=true;
        else
            leg=char(leg,'space-time');
        end
    end
    if should_I_run_this(list_runs,'brute_force_an_prec')
        plot(t_,amplitude_norm_an_prec,'+-');
        if ~has_leg
            leg=char('space-time-ANALY');
            has_leg=true;
        else
            leg=char(leg,'space-time-ANALY');
        end
    end
    if should_I_run_this(list_runs,'iqs_an_prec')
        plot(ti,amplitude_norm_iqs2,'mx-');
        plot(time_prke_iqs,power_prke_iqs,'m.-');
        if ~has_leg
            leg=char('IQS2');
            leg=char(leg,'IQS2 fine');
            has_leg=true;
        else
            leg=char(leg,'IQS2');
            leg=char(leg,'IQS2 fine');
        end
    end
    if should_I_run_this(list_runs,'iqs')
        plot(ti,amplitude_norm_iqs,'gx-');
        if ~has_leg
            leg=char('IQS');
            has_leg=true;
        else
            leg=char(leg,'IQS');
        end
    end
    if should_I_run_this(list_runs,'prke_initial_shape')
        plot(t_,Pnorm_prke,'ro-');
        if ~has_leg
            leg=char('PRKE');
            has_leg=true;
        else
            leg=char(leg,'PRKE');
        end
    end
    
    if should_I_run_this(list_runs,'prke_exact_shape')
        plot(t_,Pnorm_prkeEX,'x-');
        if ~has_leg
            leg=char('PRKE exact');
            has_leg=true;
        else
            leg=char(leg,'PRKE exact');
        end
    end
    
    if should_I_run_this(list_runs,'prke_qs_shape')
        plot(Pnorm_prkeQS,'mx-');
        if ~has_leg
            leg=char('PRKE QS');
            has_leg=true;
        else
            leg=char(leg,'PRKE QS');
        end
    end
    
    legend(leg,'Location','Best')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end