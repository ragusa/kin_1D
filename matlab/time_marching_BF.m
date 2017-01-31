function  [varargout] = time_marching_BF( dt, ntimes, u0, FUNHANDLE)

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

% initial stuff
order=1;
if strcmp(npar.method,'BDF')
    order = npar.bdf_order;
elseif strcmp(func2str(FUNHANDLE),'solve_TD_diffusion_an_prec') && strcmp(npar.an_interp_type,'lagrange')
    order=npar.interpolation_order;
end
u=zeros(length(u0),order+1); u(:,end)=u0;
t = 0:dt:ntimes*dt;
amplitude_norm=ones(ntimes+1,1);
Ptot = dat.Ptot*ones(ntimes+1,1);

dat.ode.f_beg=zeros(npar.n,1);
dat.ode.f_end=zeros(npar.n,1);


%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if io.console_print, fprintf('time end = %g \n',time_end); end
    
    u(:,1:end-1) = u(:,2:end);
    
    % solve time-dependent diffusion for flux
    if it<order;
        od = it;
    else
        od = order;
    end
    u(:,end) = FUNHANDLE(u(:,end-od:end-1),dt,t(it-od+1:it+1));
    
    % update data for hermite interpolation
    dat.ode.f_beg=dat.ode.f_end;
    
    % plot/movie
    if io.plot_transient_figure
        figure(1);
        plot(npar.x_dofs,u(1:npar.n,end));drawnow;
        if io.make_movie, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    Ptot(it+1) = compute_power(dat.nusigf,time_end,u(1:npar.n,end));
    
    % ratio of <u*,IVu> to its initial value
    amplitude_norm(it+1) = (npar.phi_adj'*npar.IV*u(1:npar.n,end))/npar.K0;
    
    % save flux for usage in PRKE exact for testing purposes
    if io.save_flux
        io.phi_save(:,it)=u(1:npar.n,end); %/Pnorm(it+1);
    end
    
    if io.print_progress
        npar.nnn = npar.nnn + 1;
        prog = npar.nnn/npar.ntot;
        progstr = sprintf('Progress: %.3g %% \n', prog*100);
        labels = {''};
        figure(101)
        pie([prog],labels);title(progstr);drawnow;
    end
end

% make movie
if io.plot_transient_figure && io.make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, movie_name, 'compression','None', 'fps',1);
end

% plot power level
if io.plot_power_figure
    plot_power_level( func2str(FUNHANDLE), ...
        linspace(0,dt*ntimes,ntimes+1), amplitude_norm);
end

if dat.PbID==12
    amplitude_norm = compute_L2norm(@(x) npar.phi_exact(x,time_end),u(1:npar.n,end))
%     amplitude_norm = abs(npar.phi_exact(npar.x_dofs(2),time_end) - u(2,end))
%     [Cend] = assemble_source(npar.C_exact,time_end);

%     amplitude_norm = norm(Cend-u(npar.n+1:end,end),2)
end
    
% output
nOutputs = nargout;
varargout = cell(1,nOutputs);
switch nOutputs
    case 1
        varargout{1} = amplitude_norm(end);
    case 2
        varargout{1} = amplitude_norm;
        varargout{2} = Ptot;
    otherwise
        error('Wrong number of output arguments in %s',mfilename);
end

end

