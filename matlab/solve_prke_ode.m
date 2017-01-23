function [X,dXdt,w,t,pow] =  solve_prke_ode(X,dt_macro,time_end,shape_beg,shape_end)

% make the problem-data a global variable
global dat npar

% create a list of times when reactivity will be recomputed for this macro
% time step
time_beg = time_end - dt_macro;

dat.ode.time_beg = time_beg;
dat.ode.time_end = time_end;

% tolerances for odesolvers
rtol = 3e-14; abso = 3e-14;
atol  = abso*ones(length(X),1);
options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);

% ODE solve for PRKEs
if npar.solve_prke_compute_rho_each_time
    % save shape to retrieve them in funprke
    dat.ode.shape_beg = shape_beg;
    dat.ode.shape_end = shape_end;
    % call ode solver with prke function that re-compute rho/beta at each time (expensive)
    [t,y]=ode15s(@funprke,[time_beg time_end],X,options);
    dXdt = funprke(time_end,y(end,:)');
    w = 0;
else
    switch npar.iqs_prke_interpolation_method
        case 1% compute prke parameters at beg/end time and linearly interpolate
            w = [ shape_beg shape_end]; 
            [dat.ode.rho_MGT_beg,dat.ode.beff_MGT_beg,dat.ode.q_MGT_beg]=compute_prke_parameters(time_beg,shape_beg);
            [dat.ode.rho_MGT_end,dat.ode.beff_MGT_end,dat.ode.q_MGT_end]=compute_prke_parameters(time_end,shape_end);
            % call ode solver with prke function that linearly interpolates
            [t,y]=ode15s(@funprke_lin_interp,[time_beg time_end],X,options);
            dXdt = funprke_lin_interp(time_end,y(end,:)');
        case 2 % linear variation for XS, linear for shape
            w = [ shape_beg shape_end]; 
            for k=1:size(w,2)
                [dat.ode.rho_MGT_beg(k),dat.ode.beff_MGT_beg(k),dat.ode.q_MGT_beg(k)]=compute_prke_parameters(time_beg,w(:,k));
                [dat.ode.rho_MGT_end(k),dat.ode.beff_MGT_end(k),dat.ode.q_MGT_end(k)]=compute_prke_parameters(time_end,w(:,k));
            end
            % call ode solver with prke function that linearly interpolates
            [t,y]=ode15s(@funprke_linlin_interp,[time_beg time_end],X,options);
            dXdt = funprke_linlin_interp(time_end,y(end,:)');
        case 3 % linear variation for XS, cubic hermite for shape
            % hermite interpolant
            t1=time_beg;
            t2=time_end;
            mat=[ t1^3 t1^2 t1 1; t2^3 t2^2 t2 1;  3*t1^2 2*t1 1 0; 3*t2^2 2*t2 1 0];
            rhs = [ shape_beg' ;shape_end' ; dat.ode.f_beg'; dat.ode.f_end'];
            % contains the w coefficients such that :
            %  shape(t) = w(1) t^3 + w(2) t^2 + w(3) t + w(4)
            w = mat\rhs; w=w'; dat.ode.nw = size(w,2);
            for k=1:dat.ode.nw
                [dat.ode.rho_MGT_beg(k),dat.ode.beff_MGT_beg(k),dat.ode.q_MGT_beg(k)]=compute_prke_parameters(time_beg,w(:,k));
                [dat.ode.rho_MGT_end(k),dat.ode.beff_MGT_end(k),dat.ode.q_MGT_end(k)]=compute_prke_parameters(time_end,w(:,k));
            end
            % call ode solver with prke function that linearly interpolates
            [t,y]=ode15s(@funprke_linhermite_interp,[time_beg time_end],X,options);
            dXdt = funprke_linhermite_interp(time_end,y(end,:)');
        case 4 % linear variation for XS, quadratic hermite for shape
            % hermite interpolant
            t1=time_beg;
            t2=time_end;
            mat=[ t1^2 t1 1; t2^2 t2 1; 2*t1 1 0];
            rhs = [ shape_beg' ;shape_end' ; dat.ode.f_beg'];
            % contains the w coefficients such that :
            %  shape(t) = w(1) t^2 + w(2) t + w(3)
            w = mat\rhs; w=w'; dat.ode.nw = size(w,2);
            for k=1:dat.ode.nw
                [dat.ode.rho_MGT_beg(k),dat.ode.beff_MGT_beg(k)]=compute_prke_parameters(time_beg,w(:,k));
                [dat.ode.rho_MGT_end(k),dat.ode.beff_MGT_end(k)]=compute_prke_parameters(time_end,w(:,k));
            end
            % call ode solver with prke function that linearly interpolates
            [t,y]=ode15s(@funprke_linhermite_interp,[time_beg time_end],X,options);
            dXdt = funprke_linhermite_interp(time_end,y(end,:)');
        otherwise
            error('interpolation_method method implemented')
    end
end
% save data
X=(y(end,:))';
pow=y(:,1);

% % save value for dpdt at the end of the macro time step
% if npar.solve_prke_compute_rho_each_time
%     dXdt=funprke(time_end,X);
% else
%     % react coef have already been computed. no need to re-compute them in
%     % the call to the prke function
%     switch npar.iqs_prke_interpolation_method
%         case 1
%             dXdt=funprke_lin_interp(time_end,X);
%         case 2
%             dXdt=funprke_linlin_interp(time_end,X);
%         case 3
%             dXdt=funprke_linhermite_interp(time_end,X);
%     end
% end
% dpdt=dXdt(1);

return
end


