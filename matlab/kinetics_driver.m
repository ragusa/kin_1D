function kinetics_driver

clc;
% clear all;
close all;
% turn off warning due to interp1
warning('OFF','MATLAB:interp1:ppGriddedInterpolant')

global npar io

% verbose/output parameters
io.console_print         = true;
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
t_end = 1.25;
ntimes=50;
dt = t_end/ntimes;
% testing with another weighting function
% npar.phi_adj = ones(length(npar.phi_adj),1);
% npar.phi_adj(1)=0;
% npar.phi_adj(end)=0;

i=0;
i=i+1; list_runs{i}= 'brute_force_matlab';
% i=i+1; list_runs{i}= 'brute_force';
% i=i+1; list_runs{i}= 'brute_force_an_prec';
i=i+1; list_runs{i}= 'brute_force_elim_prec';
% i=i+1; list_runs{i}= 'iqs_an_prec';
% i=i+1; list_runs{i}= 'iqs_an_prec1';
% i=i+1; list_runs{i}= 'iqs_elim_prec';
% i=i+1; list_runs{i}= 'iqsPC_an_prec';
i=i+1; list_runs{i}= 'iqsPC_elim_prec';
% i=i+1; list_runs{i}= 'iqs_theta_prec';
% i=i+1; list_runs{i}= 'iqs';
% npar.iqs_prke_interpolation_method=3
% i=i+1; list_runs{i}= 'prke_initial_shape';
% i=i+1; list_runs{i}= 'prke_exact_shape';
% i=i+1; list_runs{i}= 'prke_qs_shape';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB time discretization of
%%%   the TD neutron diffusion eq and precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'brute_force_matlab')
    [ref.ampl ref.Ptot] = reference_solution( t_end, u0);    
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
%%%   the TD neutron diffusion eq and precursors eq (but with elimination
%%%   of the precursors eq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'brute_force_elim_prec')
    FUNHANDLE = @solve_TD_diffusion_elim_prec;
    [brute_force_elim_prec.ampl, brute_force_elim_prec.Ptot]=time_marching_BF( dt, ntimes, u0, FUNHANDLE);
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
% if should_I_run_this(list_runs,'iqs_an_prec1')
%     npar.iqs_prke_interpolation_method=1;
%     FUNHANDLE = @solve_IQS_diffusion_an_prec;
%     [iqs_an_prec.ampl, iqs_an_prec.Ptot, iqs_an_prec.time_prke_iqs, iqs_an_prec.power_prke_iqs]=time_marching_IQS( dt, ntimes, u0, FUNHANDLE);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS IQS IQS with elimination of the precursors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqs_elim_prec')
    FUNHANDLE = @solve_IQS_diffusion_elim_prec;
    [iqs_elim_prec.ampl, iqs_elim_prec.Ptot, iqs_elim_prec.time_prke_iqs, iqs_elim_prec.power_prke_iqs]=time_marching_IQS( dt, ntimes, u0, FUNHANDLE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS version inspired from the PC-IQS method, with ANALYTICAL precursors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqsPC_an_prec')
    FUNHANDLE = @solve_IQS_PC_diffusion_an_prec;
    [iqsPC_an_prec.ampl, iqsPC_an_prec.Ptot, iqsPC_an_prec.time_prke_iqs, iqsPC_an_prec.power_prke_iqs]=time_marching_IQS( dt, ntimes, u0, FUNHANDLE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS version inspired from the PC-IQS method, with elimiation of precursors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'iqsPC_elim_prec')
    FUNHANDLE = @solve_IQS_PC_diffusion_elim_prec;
    [iqsPC_elim_prec.ampl, iqsPC_elim_prec.Ptot, iqsPC_elim_prec.time_prke_iqs, iqsPC_elim_prec.power_prke_iqs]=time_marching_IQS( dt, ntimes, u0, FUNHANDLE);
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
    FUNHANDLE = solve_PRKE;
    prke_initial_shape.ampl=time_marching_PRKE( dt, ntimes, u0, FUNHANDLE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prke using exact shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'prke_exact_shape')
    FUNHANDLE = prke_exact_shape;
    prke_exact_shape.ampl=time_marching_PRKE( dt, ntimes, u0, FUNHANDLE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prke using QS approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if should_I_run_this(list_runs,'prke_qs_shape')
    FUNHANDLE = prke_qs_shape;
    prke_qs_shape.ampl=time_marching_PRKE( dt, ntimes, u0, FUNHANDLE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
end