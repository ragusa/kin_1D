function kinetics_driver_conv

clc;
% clear all;
close all;
% turn off warning due to interp1
warning('OFF','MATLAB:interp1:ppGriddedInterpolant')

global npar io

% verbose/output parameters
io.console_print         = false;
io.plot_transient_figure = false;
io.plot_power_figure     = false;
io.make_movie            = false;
io.save_flux             = false;
io.print_progress        = false;
io.figID = 99;
% one of the two choices for applying BC
npar.set_bc_last=true;

% select problem
pbID=11; refinements=1;
problem_init(pbID,refinements);

% compute fundamental eigenmode
curr_time=0;
if pbID~=12
    [phi0]=steady_state_eigenproblem(curr_time);
else
    [phi0] = npar.phi_exact(npar.x_dofs',curr_time);
end

% initialize kinetic values
C0 = kinetics_init(phi0,curr_time);

% initial solution vector
u=[phi0;C0];
% save a copy of it
u0=u;

% time steping data
t_end = 1.28;
% t_end = 1.3;
% t_end = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB time discretization of
%%%   the TD neutron diffusion eq and precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pbID==12
    amplitude_norm_ref = 0;
else
    amplitude_norm_ref = reference_solution( t_end, u0);
end

nn=9;
ntimes = 2.^(3:nn-1)*4;
dt = t_end./ntimes;

% Interpolation type of shape for IQS prke parameters
% (see solve_prke)
%   'none' = parameters are evaluated every micro step
%   1      = linear interpolation of parameters evauated at macro steps
%   2      = linear interpolation of shape
%   3      = cubic interpolation of shape
%   4      = quadratic interpolation of shape
npar.iqs_prke_interpolation_method=2;

% Type of precursor interpolation for analytical elimination
% (see solve_TD_diffusion_an_prec)
% Redefine near where the function is called below
npar.an_interp_type='lagrange'; % This is just a place holder, actual definition is below
    
% Number of steps in history for lagrange interpolation of precursors for
% analytical elimination precursors (see solve_TD_diffusion_an_prec)
npar.interpolation_order=1;

% If this is true, precursors are recalculated using a cubic hermite
% interpolation of shape (see solve_TD_diffusion_an_prec line 192)
npar.hermite_prec_update=true;

% Used of IQS_PC, determines how precursors are solved after diffusion
% evaluation (see solve_IQS_PC_diffusion_elim_prec line 80)
% 'none'   = do not update precursor
% 'rk'     = runge-kutta revaluation
% 'linear' = linear interpolatin of shape
% 'H2'     = quadratic hermite interpolation of shape
% 'H3'     = cubic hermite interpolation of shape
npar.prec_solve_type = 'linear';


i=0;
% not to be used for conv. studies % i=i+1; list_runs{i}= 'brute_force_matlab';
% i=i+1; list_runs{i}= 'brute_force';
% i=i+1; list_runs{i}= 'brute_force_elim_prec';
i=i+1; list_runs{i}= 'brute_force_an_prec';
i=i+1; list_runs{i}= 'iqs_an_prec';
% i=i+1; list_runs{i}= 'iqs_elim_prec';
i=i+1; list_runs{i}= 'iqsPC_an_prec';
% i=i+1; list_runs{i}= 'iqsPC_elim_prec';
% i=i+1; list_runs{i}= 'iqs_theta_prec';
% i=i+1; list_runs{i}= 'iqs';
% % npar.iqs_prke_interpolation_method=3

% i=i+1; list_runs{i}= 'prke_initial_shape';
% i=i+1; list_runs{i}= 'prke_exact_shape';
% i=i+1; list_runs{i}= 'prke_qs_shape';

npar.ntot = sum(ntimes)*i;
npar.nnn = 0;

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
        display('brute_force')
        npar.an_interp_type='lagrange';
        FUNHANDLE = @solve_TD_diffusion;
        brute_force.ampl(iconv) = time_marching_BF( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% brute force discretization of
    %%%   the TD neutron diffusion eq and ANALYTICAL
    %%%   solution for the precursors eq
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'brute_force_an_prec')
        display('brute_force_an_prec')
        npar.an_interp_type='lagrange';
        FUNHANDLE = @solve_TD_diffusion_an_prec;
        brute_force_an_prec.ampl(iconv) = time_marching_BF( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% brute force discretization of
    %%%   the TD neutron diffusion eq and precursors eq (but with elimination
    %%%   of the precursors eq)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'brute_force_elim_prec')
        display('brute_force_elim_prec')
        FUNHANDLE = @solve_TD_diffusion_elim_prec;
        [brute_force_elim_prec.ampl(iconv)]=time_marching_BF( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS IQS IQS with ANALYTICAL precursors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqs_an_prec')
        display('iqs_an_prec')
%         npar.an_interp_type='hermite';
        FUNHANDLE = @solve_IQS_diffusion_an_prec;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqs_an_prec.ampl(iconv)=a;
        iqs_an_prec.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS IQS IQS with elimination of the  precursors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqs_elim_prec')
        display('iqs_elim_prec')
        FUNHANDLE = @solve_IQS_diffusion_elim_prec;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqs_elim_prec.ampl(iconv)=a;
        iqs_elim_prec.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS version inspired from the PC-IQS method, with ANALYTICAL precursors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqsPC_an_prec')
        display('iqsPC_an_prec')
        FUNHANDLE = @solve_IQS_PC_diffusion_an_prec;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqsPC_an_prec.ampl(iconv)=a;
        iqsPC_an_prec.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS version inspired from the PC-IQS method, with elimination of precursors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqsPC_elim_prec')
        display('iqsPC_elim_prec')
        FUNHANDLE = @solve_IQS_PC_diffusion_elim_prec;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqsPC_elim_prec.ampl(iconv)=a;
        iqsPC_elim_prec.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS IQS IQS with Theta Discretized precursors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqs_theta_prec')
        display('iqs_theta_prec')
        FUNHANDLE = @solve_IQS_diffusion_td_prec;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqs_theta_prec.ampl(iconv)=a;
        iqs_theta_prec.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IQS IQS IQS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'iqs')
        display('iqs')
        FUNHANDLE = @solve_IQS_diffusion;
        [a,p] = time_marching_IQS( dt(iconv), ntimes(iconv), u0, FUNHANDLE);
        iqs.ampl(iconv)=a;
        iqs.power_prke_iqs(iconv)=p;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% standard PRKE with initial shape
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'prke_initial_shape')
        display('prke_initial_shape')
        FUNNAME = solve_PRKE;
        prke_initial_shape.ampl(iconv) = time_marching_PRKE( dt(iconv), ntimes(iconv), u0, FUNNAME);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% prke using exact shape
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'prke_exact_shape')
        display('prke_exact_shape')
        FUNNAME = solve_PRKE_exact;
        prke_exact_shape.ampl(iconv) = time_marching_PRKE( dt(iconv), ntimes(iconv), u0, FUNNAME);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% prke using QS approximation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if should_I_run_this(list_runs,'prke_qs_shape')
        display('prke_qs_shape')
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
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
    end
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'brute_force_an_prec')
    error_ = abs( brute_force_an_prec.ampl - amplitude_norm_ref );
    curr_leg = 'space-time-ANALY';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
    end
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'brute_force_elim_prec')
    error_ = abs( brute_force_elim_prec.ampl - amplitude_norm_ref );
    curr_leg = 'space-time-elim';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
    end
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'iqs_an_prec')
    error_  = abs( iqs_an_prec.ampl           - amplitude_norm_ref );
    error_f = abs( iqs_an_prec.power_prke_iqs - amplitude_norm_ref );
    curr_leg = 'IQS-an';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
    end
    plot(log10(dt),log10(error_f));
    a=get(legend(gca),'String');
    leg = char(leg, char(strcat(curr_leg, ' fine')) );
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'iqs_elim_prec')
    error_  = abs( iqs_elim_prec.ampl           - amplitude_norm_ref );
    error_f = abs( iqs_elim_prec.power_prke_iqs - amplitude_norm_ref );
    curr_leg = 'IQS-elim';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
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
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
    end
    plot(log10(dt),log10(error_f));
    a=get(legend(gca),'String');
    leg = char(leg, char(strcat(curr_leg, ' fine')) );
    legend(leg,'Location','Best');
end
if should_I_run_this(list_runs,'iqsPC_elim_prec')
    error_  = abs( iqsPC_elim_prec.ampl           - amplitude_norm_ref );
    error_f = abs( iqsPC_elim_prec.power_prke_iqs - amplitude_norm_ref );
    curr_leg = 'IQS-PC-elim';
    plot(log10(dt),log10(error_),'+-');
    a=get(legend(gca),'String');
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
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
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
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
    p = polyfit(log10(dt),log10(error_),1);
    if isempty(a)
        leg=char(strcat(curr_leg,' slope = ',num2str(p(1))));
    else
        leg=char(char(a),strcat(curr_leg,' slope = ',num2str(p(1))));
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