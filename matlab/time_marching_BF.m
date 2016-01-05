function  [amplitude_norm,Ptot] = time_marching_BF( dt, ntimes, u, FUNHANDLE)

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

%%% loop on time steps %%%
for it=1:ntimes
    
    time_end=it*dt;
    if io.console_print, fprintf('time end = %g \n',time_end); end
    
    % solve time-dependent diffusion for flux
    u = FUNHANDLE(u,dt,time_end);
    
    % plot/movie
    if io.plot_transient_figure
        figure(1);
        plot(npar.x_dofs,u(1:npar.n));drawnow;
        if io.make_movie, mov(it) = getframe(gca); end
    end
    % compute end time power for plotting
    Ptot(it+1) = compute_power(dat.nusigf,time_end,u(1:npar.n));
    
    % ratio of <u*,IVu> to its initial value
    amplitude_norm(it+1) = (npar.phi_adj'*npar.IV*u(1:npar.n))/npar.K0;
    
    % save flux for usage in PRKE exact for testing purposes
%     res.phi_save(:,it)=u(1:npar.n); %/Pnorm(it+1);
    
end
% make movie
if io.plot_transient_figure && io.make_movie
    close(gcf)
    % save as AVI file
    movie2avi(mov, movie_name, 'compression','None', 'fps',1);
end

end

