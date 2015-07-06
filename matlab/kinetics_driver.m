function kinetics_driver

clear all;
close all; clc;

global dat npar

% verbose/output parameters
testing = false;
console_print = true;
plot_transient_figure = false;
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
    fprintf('%10.8g \n',keff); 
end

% initialize kinetic values
C0 = kinetics_init(phi0,curr_time);

% initial solution vector
u=[phi0;C0]; 
% save a copy of it
u0=u;

phi_adjoint = npar.phi_adj;
IV   = assemble_mass(     dat.inv_vel ,curr_time);

% time steping data
dt=0.01;
ntimes=100; % 150*2;

amplitude_norm=1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% brute force discretization of 
%%%   the TD neutron diffusion eq and precursors eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    
    TR = assemble_transient_operator(time_end);
    
    M = assemble_time_dependent_operator(time_end);
    
    % build rhs from backward Euler time discretization
    rhs = M*u;
    % build system matrix
    A = M-dt*TR;
    if npar.set_bc_last
        [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
    else
        rhs=apply_BC_vec_only(rhs);
    end
    % solve: M(unew-uold)/dt=TR.unew
    u = A\rhs;
    % plot/movie
    if plot_transient_figure
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% standard PRKE with initial shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shape=phi0;
[~,beff_MGT]=compute_prke_parameters(0.,shape);
X=[1;beff_MGT/dat.lambda];
Pnorm_prke(1)=X(1);

% loop over time steps
for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    % recompute prke parameters if XS change in time (no need for
    % up-to-date precursors, we only need rho and beta)
    [rho_MGT,beff_MGT]=compute_prke_parameters(time_end,shape);
    % build prke matrix
    A=[(rho_MGT-beff_MGT) dat.lambda ; ...
        beff_MGT         -dat.lambda];
    % solve
    X=(eye(2)-dt*A)\X;
    % store power level for plotting
    Pnorm_prke(it+1)=X(1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prke using exact shape 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shape=phi0;
[~,beff_MGT]=compute_prke_parameters(0.,shape);
X=[1;beff_MGT/dat.lambda];
Pnorm_prkeEX(1)=X(1);
IV = assemble_mass(dat.inv_vel ,0.);

% loop over time steps
for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    % compute shape from saved flux
    shape_curr = phi_save(:,it) / (phi_adjoint'*IV*phi_save(:,it)) * npar.K0;
    % recompute prke parameters  (no need for up-to-date precursors, 
    % we only need rho and beta)
    [rho_MGT,beff_MGT]=compute_prke_parameters(time_end,shape_curr);
    % build prke matrix
    A=[(rho_MGT-beff_MGT) dat.lambda ; ...
        beff_MGT         -dat.lambda];
    % solve
    X=(eye(2)-dt*A)\X;
    % store power level for plotting
    Pnorm_prkeEX(it+1)=X(1);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prke using exact QS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,beff_MGT]=compute_prke_parameters(0.,shape);
X=[1;beff_MGT/dat.lambda];
Pnorm_prkeQS(1)=X(1);
IV = assemble_mass(dat.inv_vel ,0.);

% loop over time steps
for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    % compute shape for steady state
    [shape_curr,~]=steady_state_eigenproblem(time_end);
    % recompute prke parameters  (no need for up-to-date precursors, 
    % we only need rho and beta)
    [rho_MGT,beff_MGT]=compute_prke_parameters(time_end,shape_curr);
    % build prke matrix
    A=[(rho_MGT-beff_MGT) dat.lambda ; ...
        beff_MGT         -dat.lambda];
    % solve
    X=(eye(2)-dt*A)\X;
    % store power level for plotting
    Pnorm_prkeQS(it+1)=X(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_power_figure
    figure(2); hold all;
    plot(amplitude_norm,'+-');  leg=char('space-time');
    plot(Pnorm_prkeEX,'x-');    leg=char(leg,'PRKE exact');
    plot(Pnorm_prke,'ro-');     leg=char(leg,'PRKE');
    plot(Pnorm_prkeQS,'mx-');   leg=char(leg,'PRKE QS');
    legend(leg,'Location','Best')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end