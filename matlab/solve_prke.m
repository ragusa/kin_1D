function X =  solve_prke(X,dt,time_end,shape)

% make the problem-data a global variable
global dat 

% recompute prke parameters if XS change in time 
[rho_MGT,beff_MGT]=compute_prke_parameters(time_end,shape);

% build prke matrix
A=[(rho_MGT-beff_MGT) dat.lambda ; ...
    beff_MGT         -dat.lambda];

% solve
X=(eye(2)-dt*A)\X;