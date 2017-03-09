function problem_init(problem_ID,nbr_refinements_per_region)
% load the data structure with info pertaining to the physical problem

% make the problem data a global variable
global dat npar 

% rod mov times
dat.rod_mov.t_beg_1=0.1; 
dat.rod_mov.t_end_1=0.6;
dat.rod_mov.t_beg_2=1.; 
dat.rod_mov.t_end_2=2.7;

% kinetic parameters 
dat.beta_tot=600e-5;
dat.lambda=0.1;
dat.invvel=1e-3;

% assign function pointers to the various physical coeffcients
switch problem_ID
    
    case 1
        % one material, constant in space and time   
        
        b=dat.beta_tot;
        iv=dat.invvel;
        dat.cdiff{1}   = create_material_prop('constant_in_time',1        ,[],'constant_in_space',0);
        dat.siga{1}    = create_material_prop('constant_in_time',1.0      ,[],'constant_in_space',0);
        dat.nusigf{1}  = create_material_prop('constant_in_time',1.1      ,[],'constant_in_space',0);
        dat.nusigf_p{1}= create_material_prop('constant_in_time',1.1*(1-b),[],'constant_in_space',0);
        dat.nusigf_d{1}= create_material_prop('constant_in_time',1.1*b    ,[],'constant_in_space',0);
        dat.inv_vel{1} = create_material_prop('constant_in_time',iv       ,[],'constant_in_space',0);
        dat.ext_src{1} = create_material_prop('constant_in_time',0        ,[],'constant_in_space',0);

        n_regions=1;
        region_width=400;
        dat.width = region_width * n_regions;

        imat = ones(n_regions,1);

    case 2
        % one material, constant in space      
        
        b=dat.beta_tot;
        iv=dat.invvel;
        dat.cdiff{1}   = create_material_prop('constant_in_time',1        ,[],'constant_in_space',0);
        times = [dat.rod_mov.t_beg_1 dat.rod_mov.t_end_1];
        dat.siga{1}    = create_material_prop('ramp_in_time' ,[1.0 0.98],times,'constant_in_space',0);
        dat.nusigf{1}  = create_material_prop('constant_in_time',1.1      ,[],'constant_in_space',0);
        dat.nusigf_p{1}= create_material_prop('constant_in_time',1.1*(1-b),[],'constant_in_space',0);
        dat.nusigf_d{1}= create_material_prop('constant_in_time',1.1*b    ,[],'constant_in_space',0);
        dat.inv_vel{1} = create_material_prop('constant_in_time',iv       ,[],'constant_in_space',0);
        dat.ext_src{1} = create_material_prop('constant_in_time',0        ,[],'constant_in_space',0);

        n_regions=1;
        region_width=400;
        dat.width = region_width * n_regions;

        imat = ones(n_regions,1);

    case 3
        % one material, constant in space      
        
        b=dat.beta_tot;
        iv=dat.invvel;
        dat.cdiff{1}   = create_material_prop('constant_in_time',1        ,[],'constant_in_space',0);
        times = [dat.rod_mov.t_beg_1 dat.rod_mov.t_end_1 ...
                 dat.rod_mov.t_beg_2 dat.rod_mov.t_end_2 ]/10;
        dat.siga{1}    = create_material_prop('ramp2_in_time' ,[1 0.98 1],times,'constant_in_space',0);
        dat.nusigf{1}  = create_material_prop('constant_in_time',1.1      ,[],'constant_in_space',0);
        dat.nusigf_p{1}= create_material_prop('constant_in_time',1.1*(1-b),[],'constant_in_space',0);
        dat.nusigf_d{1}= create_material_prop('constant_in_time',1.1*b    ,[],'constant_in_space',0);
        dat.inv_vel{1} = create_material_prop('constant_in_time',iv       ,[],'constant_in_space',0);
        dat.ext_src{1} = create_material_prop('constant_in_time',0        ,[],'constant_in_space',0);

        n_regions=1;
        region_width=400;
        dat.width = region_width * n_regions;

        imat = ones(n_regions,1);

    case 10
        % have material identifiers 
        n_regions = 20; % assumption: each region has the same width
        region_width=400/n_regions;
        dat.width = region_width * n_regions;
        
        imat = ones(n_regions,1);
        imat(5) = 2;
        imat(6) = 3;
        imat(15:16)=4;
        
        b=dat.beta_tot;
        iv=dat.invvel;
        dat.cdiff{1}   = create_material_prop('constant_in_time',1        ,[],'constant_in_space',0);
        dat.siga{1}    = create_material_prop('constant_in_time',1.1      ,[],'constant_in_space',0);
        dat.nusigf{1}  = create_material_prop('constant_in_time',1.1      ,[],'constant_in_space',0);
        dat.nusigf_p{1}= create_material_prop('constant_in_time',1.1*(1-b),[],'constant_in_space',0);
        dat.nusigf_d{1}= create_material_prop('constant_in_time',1.1*b    ,[],'constant_in_space',0);
        dat.inv_vel{1} = create_material_prop('constant_in_time',iv       ,[],'constant_in_space',0);
        dat.ext_src{1} = create_material_prop('constant_in_time',0        ,[],'constant_in_space',0);
        % copy material properties
        for id=2:4
            dat.cdiff{id}    = dat.cdiff{1}   ;
            dat.nusigf{id}   = dat.nusigf{1}  ;
            dat.nusigf_p{id} = dat.nusigf_p{1};
            dat.nusigf_d{id} = dat.nusigf_d{1};
            dat.inv_vel{id}  = dat.inv_vel{1} ;
            dat.ext_src{id}  = dat.ext_src{1}    ;
        end
        times = [dat.rod_mov.t_beg_1 dat.rod_mov.t_end_1];
        dat.siga{2} = create_material_prop('ramp_in_time',[1.1 1.095],times,'constant_in_space',0);
        dat.siga{4} = create_material_prop('ramp_in_time',[1.1 1.105],times,'constant_in_space',0);
        times = [dat.rod_mov.t_beg_1 dat.rod_mov.t_end_1 ...
                 dat.rod_mov.t_beg_2 dat.rod_mov.t_end_2 ];
        dat.siga{3} = create_material_prop('ramp2_in_time',[1.1 1.09 1.1],times,'constant_in_space',0);
        
    otherwise
        error('unknown problem ID ',problem_ID);
end


bc.left.type=2; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=0; % (that data is C in: -Ddu/dn=C // u/4+D/2du/dn=C // u=C)
bc.rite.type=2;
bc.rite.C=0;
dat.bc=bc; clear bc;

% normalization for nu (default=1)
npar.keff=1.;

% load the numerical parameters, npar, structure pertaining to numerics
% nbr of cells/region = nbr_refinements_per_region
% elem_to_mat = [];
% for ireg=1:dat.n_regions 
%     elem_to_mat = [elem_to_mat imat(ireg)*ones(nbr_refinements_per_region,1)];
% end
% npar.elem_to_mat = elem_to_mat; clear elem_to_mat
npar.elem_to_mat = kron(imat,ones(nbr_refinements_per_region,1));

% number of elements
npar.nel = nbr_refinements_per_region * n_regions ;
% domain
npar.x = linspace(0,dat.width,npar.nel+1);
% polynomial degree
npar.porder=1;
% nbr of dofs per variable
npar.ndofs = npar.porder*npar.nel+1;
% n: linear system size
npar.n=sum(npar.ndofs);
% ideally, we would analyze the connectivity to determine nnz
npar.nnz=(2*npar.porder+1)*npar.n; %this is an upperbound, not exact

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
poly_max=npar.porder+1;
[npar.xq,npar.wq] = GLNodeWt(poly_max);
% store shapeset
[npar.b,npar.dbdx] =feshpln(npar.xq,npar.porder);

% for newton's solve
npar.max_newton = 25;
npar.atol_newton = 1e-9;
npar.rtol_newton = 1e-9;
npar.tol_newton_lin = 1e-5;
% 1=numjac  + LU
% 2=numjac  + Gmres
% 3=numjac  + Precond Gmres
% 4=matfree + Gmres
% 5=matfree + Precond Gmres
npar.newton_solve_option = 1;
myoptP=1; 
optP=0;  if(npar.newton_solve_option==3 || npar.newton_solve_option==5), optP=myoptP; end
npar.prec_opt=optP;

% % scale variables
% npar.scale=1;
% dat.bc.left.C=dat.bc.left.C/npar.scale;
% dat.bc.rite.C=dat.bc.rite.C/npar.scale;

% connectivity
gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=[gn(iel-1,end) , gn(iel-1,2:end)+npar.porder ];
end
npar.gn=gn; clear gn;

% save x-positions of the unknowns (dofs)
npar.x_dofs=linspace(0,dat.width,npar.ndofs(1));

% save useful logical
npar.add_zero_on_diagonal=true;
npar.add_ones_on_diagonal=~npar.add_zero_on_diagonal;

return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u]=solve_ss_fem()

% perform the Newton solve
u = newton();


return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  T=newton()

global dat npar prt

% initial guess
T=ones(sum(npar.ndofs),1)*400;
T=linspace(1000,400,sum(npar.ndofs))';
T=T/npar.scale;

curr_time=0;
% compute the residual
resi = dat.myfunc(curr_time,T);
% compute its norm
norm_resi=norm(resi);
% compute stopping criterion
tol_newton = npar.atol_newton + npar.rtol_newton*norm_resi;

eigstudy=false;
if(npar.newton_solve_option>3),eigstudy=false;end

inewton=0;
while (inewton<npar.max_newton && norm_resi>tol_newton)
    inewton=inewton+1;
    if(prt.lev>=0), fprintf('Newton iteration %i, ',inewton); end
    % compute the jacobian matrix
    J = dat.myjac(curr_time,T,resi);
    P = dat.myprec(curr_time,T);
    if(eigstudy)
%         P = compute_precond(T,npar,dat,1);
        figure(2)
        subplot(2,1,1);
        eigplot(full(J),2);
        subplot(2,1,2);
        JiP=J*inv(P);
        eigplot(full(JiP),2);
        fprintf('\n\tEst. of cond number for J: %g, \t for JiP: %g\n',condest(J),condest(JiP));
        figure(3);
        subplot(2,1,1);   spy(J);
        subplot(2,1,2);   spy(JiP);
    end
    % newton solve
    switch npar.newton_solve_option
        case {1}
            dT = -J\resi;
        case{2,3}
            % gmres(A,b,restart,tol,maxit,M1,M2,x0)
            restart=[];
            tol_newton_lin=npar.tol_newton_lin;
            [dT,flag,relres,iter,resvec] = gmres(J,-resi,restart,tol_newton_lin,[],P,[],[]);
            if(prt.lev>=0 || flag==0)
                switch flag
                    case{0}
                        fprintf('\tgmres converged');
                    case{1}
                        fprintf('\tgmres did not converge');
                    case{2}
                        fprintf('\tgmres: ill-conditioned');
                    case{3}
                        fprintf('\tgmres stagnation');
                end
            end
            if(prt.lev>=0), fprintf('\titer statistics in gmres: %i, %i\n',iter(1),iter(2)); end
        case {4,5}
            restart=10;
            tol_newton_lin=npar.tol_newton_lin;
            [w,flag,relres,iter,resvec] = gmres(@Jv,-resi,restart,tol_newton_lin,150,[],[],[],... % the next line is for @afun arguments
                T,resi,dat.myfunc,P,curr_time);
            if(prt.lev>=0 || flag==0)
                switch flag
                    case{0}
                        fprintf('\tgmres converged');
                    case{1}
                        fprintf('\tgmres did not converge');
                    case{2}
                        fprintf('\tgmres: ill-conditioned');
                    case{3}
                        fprintf('\tgmres stagnation');
                end
            end
            if(prt.lev>=0), fprintf('\titer statistics in gmres: %i, %i\n',iter(1),iter(2)); end
            dT=P\w;    
    end
    % new Newton iterate
    T=T+dT;        
    % compute the residual
    resi = dat.myfunc(curr_time,T);
    % compute its norm
    norm_resi=norm(resi);
    if(prt.lev>=0)
        if(norm_resi<tol_newton)
            fprintf(' CONVERGED with ');
        end
        fprintf('Nonlinear residual error=%g, delta solution %g \n\n',norm_resi,norm(dT));
    end
end
% re-scale
T=T*npar.scale;

return
end


