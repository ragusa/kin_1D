function  [varargout] = time_marching_IQS( dt, time_final, u0, FUNHANDLE)

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

%initial stuff
if strcmp(func2str(FUNHANDLE),'solve_IQS_diffusion_an_prec') && strcmp(npar.an_interp_type,'lagrange')
    order=npar.interpolation_order;
else
    order=1;
end

u_shape=zeros(length(u0),order+1);
X = zeros(2,order+1);
tn=0;
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

% initial solution vector
 u_shape(:,end)=u0;
[~,beff_MGT]=compute_prke_parameters(0.,u0(1:npar.n));
X(:,end)=[1;beff_MGT/dat.lambda];

    dat.ode.f_beg=zeros(npar.n,1);
    dat.ode.f_end=zeros(npar.n,1);


time_prke_iqs=[];
power_prke_iqs=[];

%%% loop on time steps %%%
time_end = 0;
it = 0;
while time_end < time_final
    it = it+1;
    u_shape(:,1:end-1) = u_shape(:,2:end);
    X(:,1:end-1) = X(:,2:end);
    f_beg=dat.ode.f_beg;
    
    err = 1.0;
    while err > e_tol
        tn(it+1) = tn(it)+dt;
        time_end=tn(end);
        if io.console_print, fprintf('time end = %g \n',time_end); end
    
        % solve time-dependent diffusion for flux
        if it<order;
            od = it;
        else
            od = order;
        end

        dat.ode.f_beg=f_beg;
        
        % solve time-dependent diffusion for flux
        [u_shape(:,end),X(:,end),t,y] = FUNHANDLE(u_shape(:,end-od:end-1),X(:,end-od:end-1),dt,tn(it-od+2:it+1));
        
        if strcmp(npar.time_stepper,'DT2')
            dt_old = dt;
            [u_half,X_half,t1,y1] = FUNHANDLE(u_shape(:,end-od:end-1),X(:,end-od:end-1),dt/2,(tn(it)+tn(it+1))/2);
            dat.ode.f_beg=dat.ode.f_end;
            [u_half,X_half,t2,y2] = FUNHANDLE(u_half,X_half,dt/2,tn(it+1));
            t = [t1; t2];
            y = [y1; y2];
            err_norm = compute_L2norm(u_half(1:npar.n),u_shape(1:npar.n,end));
            shape_norm = compute_L2norm(zeros(size(u_shape(1:npar.n,end))),u_shape(1:npar.n,end));
            shape_half_norm = compute_L2norm(zeros(size(u_half(1:npar.n))),u_half(1:npar.n));
            err = err_norm/max([shape_norm shape_half_norm]) / dt;
            u_shape(:,end) = u_half;
            X(:,end) = X_half;
            dt = dt * (e_tol/err)^(1/(npar.rk.s));
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
        plot(npar.x_dofs,u_shape(1:npar.n,end));drawnow;
        if io.make_movie, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    if (strcmp(func2str(FUNHANDLE),'solve_IQS_PC_diffusion_an_prec') || ...
        strcmp(func2str(FUNHANDLE),'solve_IQS_PC_diffusion_elim_prec'))
        % u_shape is actually the flux in that case, so no need to *X(1)
        Ptot(it+1) = compute_power(dat.nusigf,time_end,         u_shape(1:npar.n,end));
    else
        Ptot(it+1) = compute_power(dat.nusigf,time_end,X(1,end)*u_shape(1:npar.n,end));
    end
    
    % ratio of <u*,IVu> to its initial value
    if (strcmp(func2str(FUNHANDLE),'solve_IQS_PC_diffusion_an_prec') || ...
        strcmp(func2str(FUNHANDLE),'solve_IQS_PC_diffusion_elim_prec'))
        % u_shape is actually the flux in that case, so no need to *X(1)
        amplitude_norm(it+1) =           (npar.phi_adj'*npar.IV*u_shape(1:npar.n,end))/npar.K0;
    else
        amplitude_norm(it+1) = X(1,end)* (npar.phi_adj'*npar.IV*u_shape(1:npar.n,end))/npar.K0;
    end
    
    % save fine-scale data
    time_prke_iqs=[time_prke_iqs ; t];
    power_prke_iqs=[power_prke_iqs; y];
    
    if io.print_progress
        npar.nnn = npar.nnn + 1;
        prog = npar.nnn/npar.ntot;
        progstr = sprintf('Progress: %.3g %% \n', prog*100);
        labels = {''};
        figure(101)
        pie([prog],labels);title(progstr); drawnow;
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
    plot_power_level( func2str(FUNHANDLE),...
        tn, amplitude_norm,...
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

