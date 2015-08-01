function X =  solve_prke(X,dt,time_end,shape)

% make the problem-data a global variable
global dat npar

% recompute prke parameters if XS change in time 
[rho_MGT,beff_MGT]=compute_prke_parameters(time_end,shape);
if npar.mf_sol
    [S_phi_MGT, S_C_MGT]=compute_mfsol_source(time_end,shape);
else
    S_phi_MGT = 0;
    S_C_MGT=0; 
end

% build prke matrix
A=[(rho_MGT-beff_MGT) dat.lambda ; ...
    beff_MGT         -dat.lambda];
rhs = X + [S_phi_MGT ; S_C_MGT]*dt;
% solve
X=(eye(2)-dt*A)\rhs;