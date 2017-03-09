function kinetics_driver

clear all;
close all; clc;

global dat npar

% verbose/output parameters
testing = false;
console_print = false;
plot_transient_figure = false;
plot_power_figure = true;
make_movie = false;

% one of the two choices for applying BC
npar.set_bc_last=true;

% select problem
pbID=10; refinements=5;
problem_init(pbID,refinements);

% compute eigenmode
curr_time=0;
[phi,keff]=steady_state_eigenproblem(curr_time);
if plot_transient_figure
    plot(npar.x_dofs,phi); 
end
if console_print
    fprintf('%10.8g \n',keff); 
end


% hold all
% [phi1,keff1]=steady_state_eigenproblem(1);
% plot(dat.x_dofs,phi1)
%
% [phi2,keff2]=steady_state_eigenproblem(3);
% plot(dat.x_dofs,phi2)
% [keff keff1 keff2]'

% initialize kinetic values
C = kinetics_init(phi,curr_time);

% initial solution vector
u=[phi;C]; 
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
phi_save=zeros(length(phi),ntimes);

% testing with another weighting function
% npar.phi_adj = ones(length(npar.phi_adj),1);
% npar.phi_adj(1)=0;
% npar.phi_adj(end)=0;

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    
    TR = assemble_transient_operator(time_end);
    
    %     load tr.mat
    %     M=A+D;
    %     M(end,end)=1;
    %     M(1,1)=1;
    %     P=NFId+NFIp;
    %     eigs(P,M,1,'lm');
    %     [uu,kk]=eigs(P,M,1,'lm');
    %     if(sum(uu)<0), uu=-uu; end
    %     flux=uu;
    %     prec=u0(npar.n+1:end);
    %     [ L*prec (NFIp-A)*flux]
    
    M = assemble_time_dependent_operator(time_end);
    
    %     load tr2.mat;
    
    % M(unew-uold)/dt=TR.unew
    rhs = M*u;
    A = M-dt*TR;
    if npar.set_bc_last
        [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
    else
        rhs=apply_BC_vec_only(rhs);
    end
    u = A\rhs;
    if plot_transient_figure, plot(npar.x_dofs,u(1:npar.n));drawnow; end
    if make_movie, mov(it) = getframe(gca); end
      
    dat.Ptot(it+1) = compute_power(dat.nusigf,time_end,u(1:npar.n));
    
    amplitude_norm(it+1) = (phi_adjoint'*IV*u(1:npar.n))/npar.K0;
    phi_save(:,it)=u(1:npar.n); %/Pnorm(it+1);
    
end
if make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, 'PbID10_v2.avi', 'compression','None', 'fps',1);
end

if testing
%%%%%%%%%%%%% test section
% end time values
u1=u;
phi1=u1(1:npar.n);
C1=u1(npar.n+1:end);

% re-load adjoint
phi_adjoint = npar.phi_adj;


% D    = assemble_stiffness(dat.cdiff   ,curr_time);
% A    = assemble_mass(     dat.siga    ,curr_time);
% NFI  = assemble_mass(     dat.nusigf  ,curr_time) / npar.keff;
% NFId = assemble_mass(     dat.nusigf_d,curr_time) / npar.keff;
% 
% IV   = apply_BC_mat_only(IV,npar.add_zero_on_diagonal);
% NFId = apply_BC_mat_only(NFId,npar.add_zero_on_diagonal);
% 
% PmM = apply_BC_mat_only(NFI-D-A,npar.add_zero_on_diagonal);
% % rho  = phi_adjoint' * (NFI-D-A) * shape;
% rho  = phi_adjoint' * (PmM) * shape;
% beff = phi_adjoint' * NFId * shape;
% MGT  = phi_adjoint' * IV * shape;

%%% check 1
% % IV   = assemble_mass(     dat.inv_vel ,curr_time);
% % 
% % shape=u0(1:npar.n);
% % aa=phi_adjoint' * IV * shape
% % norm(shape)
% % area =1/2*sum(( shape(1:end-1)+shape(2:end)));
% % aa/area
% % 
% % shape=phi1/norm(phi1)*norm(phi);
% % aa=phi_adjoint' * IV * shape 
% % norm(shape)
% % area =1/2*sum(( shape(1:end-1)+shape(2:end)));
% % aa/area
% % Pnorm=Pnorm/Pnorm(1)
% % % area =area*Pnorm(end)
% % 
%%% check 4: reduced precursors
shape=u0(1:npar.n);
[rho_MGT,beff_MGT,prec]=compute_prke_parameters(curr_time,shape,C)
IV   = assemble_mass(     dat.inv_vel ,curr_time);
MGT  = phi_adjoint' * IV * shape;
NFId = assemble_mass(     dat.nusigf_d,curr_time) / npar.keff;
bbb  = phi_adjoint' * NFId * shape / MGT
[prec (phi_adjoint'*C)/MGT]
bcc=[1 npar.n];
% shape(bcc)
% phi_adjoint(bcc)

X=[1;prec]; Xold=X;
% new shape
shape1=phi1/norm(phi_adjoint'*IV*phi1)*norm(phi_adjoint'*IV*phi);
for it=1:ntimes
    time_end=it*dt;
    [rho_MGT,beff_MGT,prec]=compute_prke_parameters(time_end,shape1,C);
    A=[(rho_MGT-beff_MGT) dat.lambda ; ...
        beff_MGT         -dat.lambda];
    X=(eye(2)-dt*A)\X
end
[X(2) (phi_adjoint'*C1)/(phi_adjoint' * IV * u0(1:npar.n))]
[X(2)/Xold(2) norm(C1)/norm(C)]
[X(2)/Xold(2) (phi_adjoint'*IV*C1)/(phi_adjoint'*IV*C)]

[   phi_adjoint' * IV * phi  ...
    phi_adjoint' * IV * shape1]
X(1)/Xold(1)
norm(phi1)/norm(phi)
norm(phi_adjoint'*IV*phi1)/norm(phi_adjoint'*IV*phi)
Pnorm/Pnorm(1)
error('end of test section');
end
%%%%%%%%%%%%% END test section

%%% standard PRKE with initial shape
shape=u0(1:npar.n);
C=u0(npar.n+1:end);
[rho_MGT,beff_MGT,prec]=compute_prke_parameters(curr_time,shape,C);
X=[1;prec];
Pnorm_prke(1)=X(1);

for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    [rho_MGT,beff_MGT,prec]=compute_prke_parameters(time_end,shape,C);
    A=[(rho_MGT-beff_MGT) dat.lambda ; ...
        beff_MGT         -dat.lambda];
    X=(eye(2)-dt*A)\X;
    Pnorm_prke(it+1)=X(1);
end


%%% prke using exact shape 
C=u0(npar.n+1:end);
[rho_MGT,beff_MGT,prec]=compute_prke_parameters(curr_time,shape,C);
X=[1;prec];
Pnorm_prkeEX(1)=X(1);
IV   = assemble_mass(     dat.inv_vel ,curr_time);

for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    shape1 = phi_save(:,it) / (phi_adjoint'*IV*phi_save(:,it)) * (phi_adjoint'*IV*phi);
    [rho_MGT,beff_MGT,prec]=compute_prke_parameters(time_end,shape1,C);
    A=[(rho_MGT-beff_MGT) dat.lambda ; ...
        beff_MGT         -dat.lambda];
    X=(eye(2)-dt*A)\X;
    Pnorm_prkeEX(it+1)=X(1);   
end

%%% prke using exact QS 
C=u0(npar.n+1:end);
[rho_MGT,beff_MGT,prec]=compute_prke_parameters(curr_time,shape,C);
X=[1;prec];
Pnorm_prkeQS(1)=X(1);
IV   = assemble_mass(     dat.inv_vel ,curr_time);

for it=1:ntimes
    time_end=it*dt;
    if console_print, fprintf('time end = %g \n',time_end); end
    [shape1,keff]=steady_state_eigenproblem(time_end);
    [rho_MGT,beff_MGT,prec]=compute_prke_parameters(time_end,shape1,C);
    A=[(rho_MGT-beff_MGT) dat.lambda ; ...
        beff_MGT         -dat.lambda];
    X=(eye(2)-dt*A)\X;
    Pnorm_prkeQS(it+1)=X(1);   
end

%%%
if plot_power_figure
    figure(2); hold all;
    plot(amplitude_norm,'+-');  leg=char('space-time');
    plot(Pnorm_prkeEX,'x-');    leg=char(leg,'PRKE exact');
    plot(Pnorm_prke,'ro-');     leg=char(leg,'PRKE');
    plot(Pnorm_prkeQS,'mx-');   leg=char(leg,'PRKE QS');
    legend(leg,'Location','Best')
end

%%%
return
end