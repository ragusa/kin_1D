function  [varargout] = time_marching_IQS( dt, ntimes, u0, FUNHANDLE)

global dat npar io

% prepare for movie
if io.make_movie
    %# figure
    figure, set(gcf, 'Color','white')
    axis([0 dat.width 0 dat.max_y_val_movie]);
    set(gca, 'nextplot','replacechildren', 'Visible','off');
    %# preallocate
    mov(1:ntimes) = struct('cdata',[], 'colormap',[]);
    % movie name
    movie_name = sprintf('PbID%d_method%s_nx%d_dt%g.avi',dat.PbID,func2str(FUNHANDLE),npar.nel,dt);
end

% initial power level
amplitude_norm=ones(ntimes+1,1);
Ptot = dat.Ptot*ones(ntimes+1,1);

% initial solution vector
u_shape=u0;
[~,beff_MGT]=compute_prke_parameters(0.,u0(1:npar.n));
X=[1;beff_MGT/dat.lambda];
if npar.iqs_prke_interpolation_method>=3
    dat.ode.f_beg=zeros(npar.n,1);
    dat.ode.f_end=zeros(npar.n,1);
end

time_prke_iqs=[];
power_prke_iqs=[];

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if io.console_print, fprintf('time end = %g \n',time_end); end
    
    % solve time-dependent diffusion for flux
    [u_shape,X,t,y] = FUNHANDLE(u_shape,X,dt,time_end);
    
    % update data for hermite interpolation
    if npar.iqs_prke_interpolation_method>=3
        dat.ode.f_beg=dat.ode.f_end;
    end
    % plot/movie
    if io.plot_transient_figure
        figure(1);
        plot(npar.x_dofs,u_shape(1:npar.n));drawnow;
        if io.make_movie, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    if (strcmp(func2str(FUNHANDLE),'solve_IQS_PC_diffusion_an_prec') || ...
        strcmp(func2str(FUNHANDLE),'solve_IQS_PC_diffusion_elim_prec'))
        % u_shape is actually the flux in that case, so no need to *X(1)
        Ptot(it+1) = compute_power(dat.nusigf,time_end,     u_shape(1:npar.n));
    else
        Ptot(it+1) = compute_power(dat.nusigf,time_end,X(1)*u_shape(1:npar.n));
    end
    
    % ratio of <u*,IVu> to its initial value
    if (strcmp(func2str(FUNHANDLE),'solve_IQS_PC_diffusion_an_prec') || ...
        strcmp(func2str(FUNHANDLE),'solve_IQS_PC_diffusion_elim_prec'))
        % u_shape is actually the flux in that case, so no need to *X(1)
        amplitude_norm(it+1) =       (npar.phi_adj'*npar.IV*u_shape(1:npar.n))/npar.K0;
    else
        amplitude_norm(it+1) = X(1)* (npar.phi_adj'*npar.IV*u_shape(1:npar.n))/npar.K0;
    end
    
    % save fine-scale data
    time_prke_iqs=[time_prke_iqs ; t];
    power_prke_iqs=[power_prke_iqs; y];
    
end
% make movie
if io.plot_transient_figure && io.make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, movie_name, 'compression','None', 'fps',1);
end

% plot power level
if io.plot_power_figure
    plot_power_level( func2str(FUNHANDLE),...
        linspace(0,dt*ntimes,ntimes+1), amplitude_norm,...
        time_prke_iqs, power_prke_iqs );
end

% output
nOutputs = nargout;
varargout = cell(1,nOutputs);
switch nOutputs
    case(2)
        varargout{1} = amplitude_norm(end);
        varargout{2} = power_prke_iqs(end);
    case(4)
        varargout{1} = amplitude_norm;
        varargout{2} = Ptot;
        varargout{3} = time_prke_iqs;
        varargout{4} = power_prke_iqs;
    otherwise
        error('Wrong number of output arguments in %s',mfilename);
end

end

