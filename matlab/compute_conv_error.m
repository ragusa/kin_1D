function [err] = compute_conv_error(phi,phiprev,fac)
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

num=0;
deno=0;

if fac < Inf
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
        fem_phi = phi(gn(iel,:));
        fem_phi_xq = b*fem_phi;
        fem_phiprev = phiprev(gn(iel,:));
        fem_phiprev_xq = b*fem_phiprev;
        num = num + Jac*dot(wq,(fem_phi_xq-fem_phiprev_xq).^fac);
        deno = deno + Jac*dot(wq,(fem_phi_xq).^fac);    
    end
    num=num^(1/fac);
    deno=deno^(1/fac);
else
    num = max(abs(phi - phiprev));
    deno = max(abs(phi));
end

err = num/deno;
return
end