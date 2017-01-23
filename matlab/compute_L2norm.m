function [err] = compute_L2norm(phiexact,phi)
% make the problem-data a global variable
global npar

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=npar.n;
% shape set
b=npar.b;
% quadrature stuff
xq=npar.xq;
wq=npar.wq;

err=0;

% loop over elements
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % jacobian of the transformation to the ref. element
    Jac=(x1-x0)/2;
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
    phiexact_xq= phiexact(x);
    fem_phi = phi(gn(iel,:));
    fem_phi_xq = b*fem_phi;
    err = err + Jac*dot(wq,(fem_phi_xq-phiexact_xq).^2);
end

err=sqrt(err);

return
end