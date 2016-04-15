function [S] = assemble_source(sfun,curr_time)

% make the problem-data a global variable
global npar

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=npar.n;
% ideally, we would analyze the connectivity to determine nnz
nnz=npar.nnz;
% allocate memory
S=zeros(n,1);
% initialize local matrices
v=zeros(porder+1,1);
% shape set
b=npar.b;
% quadrature stuff
xq=npar.xq;
wq=npar.wq;

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
    % 2/dx is the 1d jacobian
    % old: d=fct_ptr(curr_time{imat},x)/Jac;
    s = sfun(x,curr_time)*Jac;%/npar.keff;

    % assemble
    for i=1:porder+1
        v(i) = dot(s.*wq , b(:,i));
    end
    S(gn(iel,:)) = S(gn(iel,:)) + v;
end

return
end
