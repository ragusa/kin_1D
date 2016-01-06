function kinetics_driver

clc;
clear all;
close all;
% turn off warning due to interp1
warning('OFF','MATLAB:interp1:ppGriddedInterpolant')

global dat npar io res

% verbose/output parameters
io.console_print         = true;
io.plot_transient_figure = true;
io.plot_power_figure     = true;
io.make_movie            = false;

% one of the two choices for applying BC
npar.set_bc_last=true;

% select problem
pbID=10; refinements=1;
problem_init(pbID,refinements);

% compute fundamental eigenmode
curr_time=0;
[phi0]=steady_state_eigenproblem(curr_time);

% initialize kinetic values
C0 = kinetics_init(phi0,curr_time);

% initial solution vector
u=[phi0;C0];
% save a copy of it
u0=u;

% time steping data
dt=0.005;
ntimes=1.2/dt; % 150*2;
iqs_factor=1;

% save flux for tests
phi_save=zeros(length(phi0),ntimes);

% testing with another weighting function
% npar.phi_adj = ones(length(npar.phi_adj),1);
% npar.phi_adj(1)=0;
% npar.phi_adj(end)=0;

i=0;
% i=i+1; list_runs{i}= 'brute_force_matlab';
i=i+1; list_runs{i}= 'brute_force';
% i=i+1; list_runs{i}= 'brute_force_an_prec';
% i=i+1; list_runs{i}= 'iqs_an_prec';
i=i+1; list_runs{i}= 'iqsPC_an_prec';
% i=i+1; list_runs{i}= 'iqs_theta_prec';
% i=i+1; list_runs{i}= 'iqs';
% i=i+1; list_runs{i}= 'prke_initial_shape';
% i=i+1; list_runs{i}= 'prke_exact_shape';
% i=i+1; list_runs{i}= 'prke_qs_shape';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB time discretization of
%%%   the TD neutron diffusion eq and precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'brute_force_matlab')
    
    u=u0;
    
    % options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10,'OutputFcn',@odephas2,'OutputSel',[1:npar.ndofs]);
    % options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10,'OutputFcn',@odeplot,'OutputSel',[1:npar.ndofs]);
    % options = odeset('RelTol',rtol,'AbsTol',atol);
    options = odeset('RelTol',1e-10,'AbsTol',1e-10,'InitialStep',1e-8,'Stats','on','MaxStep',1e-1);
    % options = odeset('RelTol',1e-4,'AbsTol',1e-4);
    [t_ref,u_arr]=ode15s(@residual_TD_diffusion,[0 (ntimes*dt)],u0,options);
    u_arr=u_arr';
    
    %%% post-process solution %%%
    if io.plot_transient_figure
        figure;
    end
    for it=1:length(t_ref)
        % plot/movie
        if io.plot_transient_figure
            plot(npar.x_dofs,u_arr(1:npar.n,it));drawnow;
            if io.make_movie && it>1, mov(it) = getframe(gca); end
        end
        % compute end time power for plotting
        dat.Ptot(it) = compute_power(dat.nusigf,t_ref(it),u_arr(1:npar.n,it));
        % ratio of <u*,IVu> to its initial value
        amplitude_norm_ref(it) = (npar.phi_adj'*npar.IV*u_arr(1:npar.n,it))/npar.K0;
    end
    % make movie
    if io.plot_transient_figure && io.make_movie
        close(gcf)
        % save as AVI file
        movie2avi(mov, 'PbID10_v2_REF.avi', 'compression','None', 'fps',1);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% brute force discretization of
%%%   the TD neutron diffusion eq and precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'brute_force')
    FUNHANDLE = @solve_TD_diffusion;
    [brute_force.ampl, brute_force.Ptot]=time_marching_BF( dt, ntimes, u0, FUNHANDLE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% brute force discretization of
%%%   the TD neutron diffusion eq and ANALYTICAL
%%%   solution for the precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'brute_force_an_prec')
    FUNHANDLE = @solve_TD_diffusion_an_prec;
    [brute_force_an_prec.ampl, brute_force_an_prec.Ptot]=time_marching_BF( dt, ntimes, u0, FUNHANDLE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS IQS IQS with ANALYTICAL precursors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqs_an_prec')
    FUNHANDLE = @solve_IQS_diffusion_an_prec;
    [iqs_an_prec.ampl, iqs_an_prec.Ptot, iqs_an_prec.time_prke_iqs, iqs_an_prec.power_prke_iqs]=time_marching_IQS( dt, ntimes, u0, FUNHANDLE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS version inspired from the PC-IQS method, with ANALYTICAL precursors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqsPC_an_prec')
    FUNHANDLE = @solve_IQS_PC_diffusion_an_prec;
    [iqsPC_an_prec.ampl, iqsPC_an_prec.Ptot, iqsPC_an_prec.time_prke_iqs, iqsPC_an_prec.power_prke_iqs]=time_marching_IQS( dt, ntimes, u0, FUNHANDLE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS IQS IQS with Theta Discretized precursors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqs_theta_prec')
    FUNHANDLE = @solve_IQS_diffusion_td_prec;
    [iqs_theta_prec.ampl, iqs_theta_prec.Ptot, iqs_theta_prec.time_prke_iqs, iqs_theta_prec.power_prke_iqs]=time_marching_IQS( dt, ntimes, u0, FUNHANDLE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS IQS IQS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqs')
    FUNHANDLE = @solve_IQS_diffusion;
    [iqs.ampl, iqs.Ptot, iqs.time_prke_iqs, iqs.power_prke_iqs]=time_marching_IQS( dt, ntimes, u0, FUNHANDLE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% standard PRKE with initial shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'prke_initial_shape')
    
    shape0=phi0;
    [~,beff_MGT]=compute_prke_parameters(0.,shape0);
    X=[1;beff_MGT/dat.lambda];
    Pnorm_prke=X(1)*ones(ntimes+1);
    
    % loop over time steps
    for it=1:ntimes
        time_end=it*dt;
        if io.console_print, fprintf('time end = %g \n',time_end); end
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
    Pnorm_prkeEX=X(1)*ones(ntimes+1);
    IV = assemble_mass(dat.inv_vel ,0.);
    
    % loop over time steps
    for it=1:ntimes
        time_end=it*dt;
        if io.console_print, fprintf('time end = %g \n',time_end); end
        % compute shape from saved flux
        shape_curr = phi_save(:,it) / (npar.phi_adj'*npar.IV*phi_save(:,it)) * npar.K0;
        % solve prke
        X =  solve_prke(X,dt,time_end,shape_curr);
        % store power level for plotting
        Pnorm_prkeEX(it+1)=X(1);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prke using QS approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'prke_qs_shape')
    
    shape0=phi0;
    [~,beff_MGT]=compute_prke_parameters(0.,shape0);
    X=[1;beff_MGT/dat.lambda];
    Pnorm_prkeQS=X(1)*ones(ntimes+1);
    IV = assemble_mass(dat.inv_vel ,0.);
    
    % loop over time steps
    for it=1:ntimes
        time_end=it*dt;
        if io.console_print, fprintf('time end = %g \n',time_end); end
        % compute shape for steady state
        [shape_curr,~]=steady_state_eigenproblem(time_end);
        % solve prke
        X =  solve_prke(X,dt,time_end,shape_curr);
        % store power level for plotting
        Pnorm_prkeQS(it+1)=X(1);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if io.plot_power_figure
    figure(3); hold all;
    t_=linspace(0,dt*ntimes,ntimes+1);
    ti=linspace(0,dt*ntimes,ntimes/iqs_factor+1);
    has_leg=false;
    
    if should_I_run_this(list_runs,'brute_force_matlab')
        plot(t_ref,amplitude_norm_ref,'+-');
        if ~has_leg
            leg=char('space-time-matlab');
            has_leg=true;
        else
            leg=char(leg,'space-time-matlab');
        end
    end
    if should_I_run_this(list_runs,'brute_force')
        plot(t_,brute_force.amplitude_norm,'+-');
        if ~has_leg
            leg=char('space-time');
            has_leg=true;
        else
            leg=char(leg,'space-time');
        end
    end
    if should_I_run_this(list_runs,'brute_force_an_prec')
        plot(t_,brute_force_an_prec.amplitude_norm,'+-');
        if ~has_leg
            leg=char('space-time-ANALY');
            has_leg=true;
        else
            leg=char(leg,'space-time-ANALY');
        end
    end
    if should_I_run_this(list_runs,'iqs_an_prec')
        plot(ti,iqs_an_prec.amplitude_norm,'mx-');
        plot(time_prke_iqs,power_prke_iqs,'m.-');
        if ~has_leg
            leg=char('IQS-an');
            leg=char(leg,'IQS-an fine');
            has_leg=true;
        else
            leg=char(leg,'IQS-an');
            leg=char(leg,'IQS-an fine');
        end
    end
    if should_I_run_this(list_runs,'iqsPC_an_prec')
        plot(ti,iqsPC_an_prec.amplitude_norm,'o-');
        plot(time_prke_iqs_pc,power_prke_iqs_pc,'c.-');
        if ~has_leg
            leg=char('IQS-PC-an');
            leg=char(leg,'IQS-PC-an fine');
            has_leg=true;
        else
            leg=char(leg,'IQS-PC-an');
            leg=char(leg,'IQS-PC-an fine');
        end
    end
    if should_I_run_this(list_runs,'iqs_theta_prec')
        plot(ti,iqs_theta_prec.amplitude_norm,'kx-');
        if ~has_leg
            leg=char('IQS-theta-prec');
            has_leg=true;
        else
            leg=char(leg,'IQS-theta-prec');
        end
    end
    if should_I_run_this(list_runs,'iqs')
        plot(ti,iqs.amplitude_norm,'gx-');
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