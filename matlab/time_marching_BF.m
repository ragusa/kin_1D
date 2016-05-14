function  [varargout] = time_marching_BF( dt, time_final, u0, FUNHANDLE)

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
if strcmp(func2str(FUNHANDLE),'solve_TD_diffusion_an_prec') && strcmp(npar.an_interp_type,'lagrange') && strcmp(npar.time_stepper,'constant')
    order=npar.interpolation_order;
else
    order=1;
end

u=zeros(length(u0),order+1); u(:,end)=u0;
t=0;
if strcmp(npar.time_stepper,'constant')
    ntimes = time_final/dt;
%     t = 0:dt:ntimes*dt;
    amplitude_norm=ones(ntimes+1,1);
    Ptot = dat.Ptot*ones(ntimes+1,1);
    e_tol = 0.1;
else
    amplitude_norm=ones(2,1);
    Ptot = dat.Ptot*ones(2,1);
    e_tol = npar.time_tol;
end

dat.ode.f_beg=zeros(npar.n,1);
dat.ode.f_end=zeros(npar.n,1);


%%% loop on time steps %%%
time_end = 0;
it = 0;
while time_end < time_final
    it = it+1;
    u(:,1:end-1) = u(:,2:end);
    f_beg=dat.ode.f_beg;
    
    err = 1.0;
    while err > e_tol
        t(it+1) = t(it)+dt;
        time_end=t(end);
        if io.console_print, fprintf('time end = %g \n',time_end); end    

        % solve time-dependent diffusion for flux
        if it<order;
            od = it;
        else
            od = order;
        end
                
        dat.ode.f_beg=f_beg;
        
        u(:,end) = FUNHANDLE(u(:,end-od:end-1),dt,t(it-od+2:it+1));

        if strcmp(npar.time_stepper,'DT2')
            dt_old = dt;
            u_half = FUNHANDLE(u(:,end-od:end-1),dt/2,(t(it)+t(it+1))/2);
            dat.ode.f_beg=dat.ode.f_end;
            u_half = FUNHANDLE(u_half,dt/2,t(it+1));
            err_norm = compute_L2norm(u_half(1:npar.n),u(1:npar.n,end));
            phi_norm = compute_L2norm(zeros(size(u(1:npar.n,end))),u(1:npar.n,end));
            phi_half_norm = compute_L2norm(zeros(size(u_half(1:npar.n))),u_half(1:npar.n));
            err = err_norm/max([phi_norm phi_half_norm]);
            u(:,end) = u_half;
            dt = dt * (e_tol/err)^(1/npar.rk.s);
            if err <= e_tol
                if time_final<time_end+dt
                    dt = time_final - time_end;
                end
                if dt/dt_old > npar.max_increase
                    dt = npar.max_increase*dt_old;
                end
            end
            
        else
            err = 0.0;
        end
    end
        
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
    plot_power_level( func2str(FUNHANDLE), t, amplitude_norm);
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

