function [S] = assemble_source(sfun,curr_time)

porder= 3;
gn    = [1,2,3,4];
n=4;
% allocate memory
S=zeros(n,1);
% initialize local matrices
v=zeros(porder+1,1);
% shape set
b=[0.660005665072803,0.520937687711704,-0.230187903250739,0.0492445504662318;...
   0.00337373643277257,1.00488585482565,-0.00992135357232475,0.00166176231390641;...
   0.00166176231390640,-0.00992135357232467,1.00488585482565,0.00337373643277254;...
   0.0492445504662318,-0.230187903250739,0.520937687711703,0.660005665072804];
% quadrature stuff
xq=[-0.861136311594053;-0.339981043584856;0.339981043584856;0.861136311594053];
wq=[0.347854845137454;0.652145154862546;0.652145154862546;0.347854845137454];

% element extremities
x0=0;
x1=1;
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
S(gn(1,:)) = S(gn(1,:)) + v;

