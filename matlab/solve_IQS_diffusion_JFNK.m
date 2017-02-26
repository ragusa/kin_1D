function [u_shape, X,t,y] = solve_IQS_diffusion_JFNK(u_shape,X,dt_macro,tn)

global io dat npar

% shortcuts
dt = dt_macro;
C_old = u_shape(npar.n+1:end,end);

% save values at beginning of macro time step: they are needed in the IQS iteration
X_beg=X;
% for clarity, I also single out the shape function at the beginning/end of the macro time step
shape_beg=u_shape(1:npar.n,:);
shape_end=shape_beg(:,end);

time_end = tn(end);
time_beg = time_end-dt;

% assemble matrices
D    = assemble_stiffness(dat.cdiff   ,time_end);
A    = assemble_mass(     dat.siga    ,time_end);
NFIp = assemble_mass(     dat.nusigf_p,time_end) / npar.keff;
IV   = assemble_mass(     dat.inv_vel ,time_end);
NFId_new = assemble_mass(dat.nusigf_d ,time_end) / npar.keff;
NFId_old = assemble_mass(dat.nusigf_d ,time_beg) / npar.keff;
TR=NFIp -(D+A);

eps = 1e-6;

for iter = 1: npar.max_iter_iqs   

    % Compute residual and precursors
    dat.solve_prec = true;
    Re = residual_IQS(shape_end,shape_beg,X_beg,tn,TR,NFId_new,NFId_old,IV,C_old);
    residual = norm(Re,2);
    if io.console_print
        fprintf('  IQS iter %d, Nonlinear|R| %g \n',iter,residual);
    end
    if residual < npar.tol_iqs
        break;
    end
    
    % Solve linear system
    dat.solve_prec = false;
    dshape = gmres(@(y) IQSfun(y,shape_end,shape_beg,X_beg,tn,TR,NFId_new,NFId_old,IV,Re,eps,C_old), -Re,[],1e-4,npar.n);
    shape_end = shape_end + dshape;
    
end

u_shape = [ shape_end ; dat.C_new];
t = dat.t;
y = dat.y;
X = dat.X;

% dat.iter_IQS = [dat.iter_IQS iter];
% dat.iter_err = [dat.iter_err err];

% renormalize anyway
if npar.scale_IQS
    u_shape = u_shape / ( ((npar.phi_adj)'*npar.IV*shape_end) / npar.K0 );
end