function S = assemble_source(fct_S,curr_time,phi)

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
a=zeros(1,porder+1);
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
    % old: s=fct_ptr(curr_time,x)*Jac;
    imat=npar.elem_to_mat(iel);
    s = fct_S(imat,x,curr_time,phi)*Jac;
    
    S(gn(iel,:)) = S(gn(iel,:)) + b(:,1).*s.*wq;
end

return
end