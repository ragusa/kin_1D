function problem_init(problem_ID,nbr_refinements_per_region)
% load the data structure with info pertaining to the physical problem

% make the problem data a global variable
global dat npar

% save pb ID for later use
dat.PbID = problem_ID;

% rod mov times
dat.rod_mov.t_beg_1=0.1; dat.rod_mov.t_end_1=0.6;
dat.rod_mov.t_beg_2=1.0; dat.rod_mov.t_end_2=2.7;

% kinetic parameters
dat.beta_tot=600e-5;
dat.lambda=0.1;
dat.invvel=1e-3;

dat.source_phi = @(x,t) 0;

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
        % one material, constant in space, ramp in time 
        
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
        % one material, constant in space, ramp2 in time 
        
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
        % have material identifiers: 20 regions. 4 materials. 3 rod movements
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
        % copy material properties that remain unchanged
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
        
    case 11
        % have material identifiers: 20 regions. 2 materials. 1 rod movement
        n_regions = 20; % assumption: each region has the same width
        region_width=400/n_regions;
        dat.width = region_width * n_regions;
        
        imat = ones(n_regions,1);
        imat(5) = 2;
        
        b=dat.beta_tot;
        iv=dat.invvel;
        dat.cdiff{1}   = create_material_prop('constant_in_time',1        ,[],'constant_in_space',0);
        dat.siga{1}    = create_material_prop('constant_in_time',1.1      ,[],'constant_in_space',0);
        dat.nusigf{1}  = create_material_prop('constant_in_time',1.1      ,[],'constant_in_space',0);
        dat.nusigf_p{1}= create_material_prop('constant_in_time',1.1*(1-b),[],'constant_in_space',0);
        dat.nusigf_d{1}= create_material_prop('constant_in_time',1.1*b    ,[],'constant_in_space',0);
        dat.inv_vel{1} = create_material_prop('constant_in_time',iv       ,[],'constant_in_space',0);
        dat.ext_src{1} = create_material_prop('constant_in_time',0        ,[],'constant_in_space',0);
        % copy material properties that remain unchanged
        for id=2:max(imat) 
            dat.cdiff{id}    = dat.cdiff{1}   ;
            dat.nusigf{id}   = dat.nusigf{1}  ;
            dat.nusigf_p{id} = dat.nusigf_p{1};
            dat.nusigf_d{id} = dat.nusigf_d{1};
            dat.inv_vel{id}  = dat.inv_vel{1} ;
            dat.ext_src{id}  = dat.ext_src{1}    ;
        end
        times = [dat.rod_mov.t_beg_1 dat.rod_mov.t_beg_2];
        dat.siga{2} = create_material_prop('ramp_in_time',[1.1 1.08],times,'constant_in_space',0);
        
    case 12
        % manufactured solution, one material, constant in space and time,
        % no precursors
        
        b=dat.beta_tot*0;
        iv=dat.invvel;
        cdiff = 1.0;
        Sa = 1.0;
        nfSf = 1.1;        
        dat.cdiff{1}   = create_material_prop('constant_in_time',cdiff    ,[],'constant_in_space',0);
        dat.siga{1}    = create_material_prop('constant_in_time',Sa       ,[],'constant_in_space',0);
        dat.nusigf{1}  = create_material_prop('constant_in_time',nfSf     ,[],'constant_in_space',0);
        dat.nusigf_p{1}= create_material_prop('constant_in_time',1.1*(1-b),[],'constant_in_space',0);
        dat.nusigf_d{1}= create_material_prop('constant_in_time',1.1*b    ,[],'constant_in_space',0);
        dat.inv_vel{1} = create_material_prop('constant_in_time',iv       ,[],'constant_in_space',0);
        dat.ext_src{1} = create_material_prop('constant_in_time',0        ,[],'constant_in_space',0);
        dat.mass{1}   = create_material_prop('constant_in_time',1.0       ,[],'constant_in_space',0);
        
        n_regions=1;
        region_width=1;
        dat.width = region_width * n_regions;
        
        imat = ones(n_regions,1);
        
        syms x t
        syms S_p(x,t) phi(x,t)
        
        syms f(t) a(x,t) A(t)
%         f(t) = (t+1)^5
%         a(x,t) = sin(x*(1-x)*(t+1))
%         A(t) = int(a(x,t),x,0,1)
%         phi(x,t) = f(t)*a(x,t)/A(t);
        f(t) = (t+1)^2;
        a(x,t) = (1+t)*x*(1-x)*(x+t);
        A(t) = int(a(x,t),x,0,1);
        phi(x,t) = f(t)*a(x,t)/A(t);
        
%         phi(x,t) = x*(1-x)*(1+t)^2;
        npar.phi_exact = matlabFunction(phi);
        
        S_p(x,t) = iv*diff(phi,t) + (Sa-nfSf)*phi - cdiff*diff(diff(phi,x),x)
        dat.source_phi = matlabFunction(S_p);
        clear x t S_p phi
        
        
    otherwise
        error('unknown problem ID ',problem_ID);
end


bc.left.type=2; % 0=neumann, 1=robin, 2=dirichlet
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

% IQS options
npar.solve_prke_compute_rho_each_time = false;
npar.prke_solve = 'matlab' ;
npar.int_order = 3;
% npar.prke_solve = 'no' ;
npar.max_iter_iqs = 6;
npar.tol_iqs      = 1e-11;
npar.iqs_prke_interpolation_method=2;

if ~strcmpi(npar.prke_solve,'matlab')
    npar.n_micro=10;
    npar.freq_react=1;
end

% movie options
dat.max_y_val_movie = 2.;

% time integration
time_integration=3;
switch time_integration
    case 1
        nstages=1;
        a=zeros(nstages,nstages);
        b=zeros(nstages,1); c=b;
        a(1,1)=1;
        c=sum(a,2);
        b(1)=1;
        
    case 3
        nstages=3;
        g=0.43586652150845899941601945119356;
        a=[g 0 0; ...
            ((1-g)/2) g 0;...
            (-(6*g^2-16*g+1)/4) ((6*g^2-20*g+5)/4) g];
        c=sum(a,2);
        b=a(nstages,:);
    otherwise
        error('time integration %d not yet implemented',time_integration);
end
npar.rk.s=nstages;
npar.rk.a=a; npar.rk.b=b; npar.rk.c=c;

return
end






