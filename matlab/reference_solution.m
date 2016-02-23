function  [varargout] = reference_solution( t_end, u0)

global dat npar io

% options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10,'OutputFcn',@odephas2,'OutputSel',[1:npar.ndofs]);
% options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10,'OutputFcn',@odeplot,'OutputSel',[1:npar.ndofs]);
% options = odeset('RelTol',rtol,'AbsTol',atol);
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'InitialStep',1e-8,'Stats','on','MaxStep',1e-1,'Stats','off');
% options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t_ref,u_arr]=ode15s(@residual_TD_diffusion,[0 t_end],u0,options);
u_arr=u_arr';


% prepare for movie
if io.make_movie
    %# figure
    figure, set(gcf, 'Color','white')
    axis([0 dat.width 0 dat.max_y_val_movie]);
    set(gca, 'nextplot','replacechildren', 'Visible','off');
    %# preallocate
    mov(1:ntimes) = struct('cdata',[], 'colormap',[]);
    % movie name
    movie_name = sprintf('PbID%d_method%s_nx%d.avi',dat.PbID,'ref_matlab',npar.nel);
end

% initial power level
amplitude_norm=ones(length(t_ref),1);
Ptot = dat.Ptot*ones(length(t_ref),1);

%%% post-process solution %%%
if io.plot_transient_figure
    figure;
end

for it=1:length(t_ref)
    % plot/movie
    if io.plot_transient_figure
        plot(npar.x_dofs,u_arr(1:npar.n,it));drawnow;
        if io.make_movie && it>1, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    Ptot(it) = compute_power(dat.nusigf,t_ref(it),u_arr(1:npar.n,it));
    % ratio of <u*,IVu> to its initial value
    amplitude_norm(it) = (npar.phi_adj'*npar.IV*u_arr(1:npar.n,it))/npar.K0;
end

% make movie
if io.plot_transient_figure && io.make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, movie_name, 'compression','None', 'fps',1);
end

% plot power level
if io.plot_power_figure
    plot_power_level( 'brute_force_matlab', t_ref, amplitude_norm);
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

