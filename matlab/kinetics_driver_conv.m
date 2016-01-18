function kinetics_driver_conv

clc;
clear all;
close all;
% turn off warning due to interp1
warning('OFF','MATLAB:interp1:ppGriddedInterpolant')

global npar io

% verbose/output parameters
io.console_print         = false;
io.plot_transient_figure = true;
io.plot_power_figure     = true;
io.make_movie            = false;
io.save_flux             = false;
io.figID = 99;
% one of the two choices for applying BC
npar.set_bc_last=true;

% select problem
pbID=11; refinements=1;
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
t_end = 1.28;
t_end = 1.3;
t_end = 1.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB time discretization of
%%%   the TD neutron diffusion eq and precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amplitude_norm_ref = reference_solution( t_end, u0);

nn=3;
ntimes = 25 * 2.^(linspace(0,nn,nn+1));
dt = t_end./ntimes;

i=0;
% % i=i+1; list_runs{i}= 'brute_force_matlab';
i=i+1; list_runs{i}= 'brute_force';
i=i+1; list_runs{i}= 'brute_force_elim_prec';
i=i+1; list_runs{i}= 'brute_force_an_prec';
% i=i+1; list_runs{i}= 'iqs_an_prec';
% i=i+1; list_runs{i}= 'iqsPC_an_prec';
% i=i+1; list_runs{i}= 'iqs_theta_prec';
% i=i+1; list_runs{i}= 'iqs';
% i=i+1; list_runs{i}= 'prke_initial_shape';
% i=i+1; list_runs{i}= 'prke_exact_shape';
% i=i+1; list_runs{i}= 'prke_qs_shape';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% LOOP OVER DT's
for iconv=1:length(ntimes)
    
    fprintf('%d out of %d\n',iconv,length(ntimes))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% brute force discretization of
    %%%   the TD neutron diffusion eq and precursors eq
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'brute_force')
        FUNHANDLE = @solve_TD_diffusion;
        brute_force.ampl(iconv) = time_marching_BF( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% brute force discretization of
    %%%   the TD neutron diffusion eq and ANALYTICAL
    %%%   solution for the precursors eq
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'brute_force_an_prec')
        FUNHANDLE = @solve_TD_diffusion_an_prec;
        brute_force_an_prec.ampl(iconv) = time_marching_BF( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% brute force discretization of
    %%%   the TD neutron diffusion eq and precursors eq (but with elimination
    %%%   of the precursors eq)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'brute_force_elim_prec')
        FUNHANDLE = @solve_TD_diffusion_elim_prec;
        [brute_force_elim_prec.ampl(iconv)]=time_marching_BF( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS IQS IQS with ANALYTICAL precursors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqs_an_prec')
        FUNHANDLE = @solve_IQS_diffusion_an_prec;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqs_an_prec.ampl(iconv)=a;
        iqs_an_prec.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS version inspired from the PC-IQS method, with ANALYTICAL precursors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqsPC_an_prec')
        FUNHANDLE = @solve_IQS_PC_diffusion_an_prec;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqsPC_an_prec.ampl(iconv)=a;
        iqsPC_an_prec.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS IQS IQS with Theta Discretized precursors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqs_theta_prec')
        FUNHANDLE = @solve_IQS_diffusion_td_prec;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqs_theta_prec.ampl(iconv)=a;
        iqs_theta_prec.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS IQS IQS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqs')
        FUNHANDLE = @solve_IQS_diffusion;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqs.ampl(iconv)=a;
        iqs.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% standard PRKE with initial shape
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'prke_initial_shape')
        FUNNAME = solve_PRKE;
        prke_initial_shape.ampl(iconv) = time_marching_PRKE( dt(iconv), ntimes(iconv), u0, FUNNAME);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% prke using exact shape
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'prke_exact_shape')
        FUNNAME = solve_PRKE_exact;
        prke_exact_shape.ampl(iconv) = time_marching_PRKE( dt(iconv), ntimes(iconv), u0, FUNNAME);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% prke using QS approximation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'prke_qs_shape')
        FUNNAME = solve_PRKE_QS;
        prke_qs_shape.ampl(iconv) = time_marching_PRKE( dt(iconv), ntimes(iconv), u0, FUNNAME);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(100); hold all;
if should_I_run_this(list_runs,'brute_force')
    error_ = abs( brute_force.ampl - amplitude_norm_ref );
    curr_leg = 'space-time';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    if isempty(a)
        leg=char(curr_leg);
    else
        leg=char(char(a),curr_leg);
    end
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'brute_force_an_prec')
    error_ = abs( brute_force_an_prec.ampl - amplitude_norm_ref );
    curr_leg = 'space-time-ANALY';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    if isempty(a)
        leg=char(curr_leg);
    else
        leg=char(char(a),curr_leg);
    end
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'brute_force_elim_prec')
    error_ = abs( brute_force_elim_prec.ampl - amplitude_norm_ref );
    curr_leg = 'space-time-elim';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    if isempty(a)
        leg=char(curr_leg);
    else
        leg=char(char(a),curr_leg);
    end
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'iqs_an_prec')
    error_  = abs( iqs_an_prec.ampl           - amplitude_norm_ref );
    error_f = abs( iqs_an_prec.power_prke_iqs - amplitude_norm_ref );
    curr_leg = 'IQS-an';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    if isempty(a)
        leg=char(curr_leg);
    else
        leg=char(char(a),curr_leg);
    end
    plot(log10(dt),log10(error_f));
    a=get(legend(gca),'String');
    leg = char(leg, char(strcat(curr_leg, ' fine')) );
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'iqsPC_an_prec')
    error_  = abs( iqsPC_an_prec.ampl           - amplitude_norm_ref );
    error_f = abs( iqsPC_an_prec.power_prke_iqs - amplitude_norm_ref );
    curr_leg = 'IQS-PC-an';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    if isempty(a)
        leg=char(curr_leg);
    else
        leg=char(char(a),curr_leg);
    end
    plot(log10(dt),log10(error_f));
    a=get(legend(gca),'String');
    leg = char(leg, char(strcat(curr_leg, ' fine')) );
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'iqs_theta_prec')
    error_  = abs( iqs_theta_prec.ampl           - amplitude_norm_ref );
    error_f = abs( iqs_theta_prec.power_prke_iqs - amplitude_norm_ref );
    curr_leg = 'IQS-theta-prec';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    if isempty(a)
        leg=char(curr_leg);
    else
        leg=char(char(a),curr_leg);
    end
    plot(log10(dt),log10(error_f));
    a=get(legend(gca),'String');
    leg = char(leg, char(strcat(curr_leg, ' fine')) );
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'iqs')
    error_  = abs( iqs.ampl           - amplitude_norm_ref );
    error_f = abs( iqs.power_prke_iqs - amplitude_norm_ref );
    curr_leg = 'IQS';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    if isempty(a)
        leg=char(curr_leg);
    else
        leg=char(char(a),curr_leg);
    end
    plot(log10(dt),log10(error_f));
    a=get(legend(gca),'String');
    leg = char(leg, char(strcat(curr_leg, ' fine')) );
    legend(leg,'Location','Best');
end
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
end