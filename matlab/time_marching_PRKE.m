function  [Pnorm_prke] = time_marching_PRKE( dt, ntimes, u, FUNHANDLE)

global dat npar io

shape0=u(npar.n);
[~,beff_MGT]=compute_prke_parameters(0.,shape0);
X=[1;beff_MGT/dat.lambda];
Pnorm_prke=X(1)*ones(ntimes+1);

% loop over time steps
for it=1:ntimes
    time_end=it*dt;
    if io.console_print, fprintf('time end = %g \n',time_end); end
    switch FUNHANDLE
        case strcmp('solve_PRKE',FUNHANDLE)
            shape_curr = shape0;
            
        case strcmp('solve_PRKE_exact',FUNHANDLE)
            if io.save_flux
                % shape from saved flux
                shape_curr = io.phi_save(:,it+1) / (npar.phi_adj'*npar.IV*io.phi_save(:,it+1)) * npar.K0;
            else
                warning('Cannot retrieve saved flux, using initial shape instead');
                shape_curr = shape0;
            end
            
        case strcmp('solve_PRKE_QS',FUNHANDLE)
            % compute shape for steady state
            shape_curr = steady_state_eigenproblem(time_end);
    end
    % solve prke
    X =  solve_prke(X,dt,time_end,shashape_currpe0);
    % store power level for plotting
    Pnorm_prke(it+1)=X(1);
end

if io.plot_power_figure
    plot_power_level( func2str(FUNHANDLE), ...
        linspace(0,dt*ntimes,ntimes+1), Pnorm_prke);
end

end

