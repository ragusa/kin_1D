function [amplitude_norm_iqs] = IQS_driver(phi0,C0,dt,ntimes,iqs_factor, plot_transient_figure,theta_log,fig)

global npar dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IQS IQS IQS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amplitude_norm_iqs=1;
theta=0;

% initial solution vector
u_shape=[phi0;C0];
[~,beff_MGT]=compute_prke_parameters(0.,phi0);
X=[1;beff_MGT/dat.lambda];
Pnorm_prkeIQS(1)=X(1);

dt=dt*iqs_factor; ntimes=ntimes/iqs_factor;

n_micro=10;
freq_react=1;

%%% loop on time steps %%%
% dat.dlogp_dt = 0;
% Lambda = 0;
% rho_minus_beta = 0;
% beta = 0;
% Lambda_orig = 0;
% time = 0;
for it=1:ntimes
    
    time_end=it*dt;
    u_shape_old=u_shape;
%     if console_print, 
        fprintf('time end = %g \n',time_end); %end
    
    % solve time-dependent diffusion for flux
    [u_shape,X,IV,err] = solve_IQS_diffusion(u_shape,X,dt,time_end,n_micro,freq_react,theta,theta_log);
    
    % plot/movie
    if plot_transient_figure
        % loop on micro step
        % inside loop, linearly interpol shape
        figure(fig)
%         subplot(2,1,1)
        if npar.mf_sol
            plot(npar.x_dofs,u_shape(1:npar.n),'b-',npar.x_dofs,1.5*npar.x_dofs/dat.width.*(1-npar.x_dofs/dat.width),'k+');drawnow;
        else
            plot(npar.x_dofs,u_shape(1:npar.n)); drawnow;
        end
%         subplot(2,1,2)
%         hold on
%         plot(time_end,err,'b+')
%         hold off
    end
    % compute end time power for plotting
    dat.Ptot_iqs(it+1) = compute_power(dat.nusigf,time_end,X(1)*u_shape(1:npar.n));
    
    % ratio of <u*,IVu> to its initial value
    amplitude_norm_iqs(it+1) = X(1)* (npar.phi_adj'*npar.IVel*u_shape(1:npar.n))/npar.K0;

%     Lambda = [Lambda; dat.Lambda_micro];
%     rho_minus_beta = [rho_minus_beta; dat.rho_minus_beta_micro];
%     beta = [beta; dat.beta_micro];
%     time = [time; time_end];
end
% Lambda_orig = ones(size(time))*npar.K0;
% dlogp_dt= dat.dlogp_dt;
% T = table(time, beta, Lambda, Lambda_orig, rho_minus_beta, dlogp_dt)