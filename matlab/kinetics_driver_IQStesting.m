function kinetics_driver_IQStesting

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
pbID=2; refinements=6;
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
% t_end = 1.28;
t_end = 1.3;
% t_end = 0.1;
ntimes=52;
dt = t_end/ntimes;
t = 0:dt:t_end;
        
% Solving for flux using BF, only capable for rk precursor elimination method
FUNHANDLE = @solve_TD_diffusion_elim_prec;
[brute_force_elim_prec.ampl,~,u_bf]=time_marching_BF( dt, ntimes, u0, FUNHANDLE);

% New function that uses the flux from the brute force method and outputs
% the amplitude.
[p_IQS] = IQS_testing(dt,u_bf);

[p_an] = analytical_prke(u_bf,dt,ntimes);

% Compare values side by side
show = [t', p_IQS, p_an]
    
figure(1)
plot(t,p_IQS,t,p_an)
legend('IQS','Analytical')
xlabel('time'); ylabel('p');

% Error between methods
error = [t', abs(p_IQS-p_an)./p_an]
    
end
