function Re = residual_IQS(shape,shape_beg,X_beg,tn,TR,NFId_new,NFId_old,IV,C_old)

global dat npar

lambda = dat.lambda;

order = length(tn)-1;
dt = tn(end) - tn(end-1);
time_end = tn(end);
time_beg = time_end-dt;

[X,dXdt,~,t,y] =  solve_prke_ode(X_beg,dt,tn,shape_beg,shape);

% interpolating polynomial
pp = interp1(t,y,'spline','pp');

p = X(1);
dpdt = dXdt(1);

% build rhs
a = zeros(order+1,2);
for j=1:(order+1)
    Aj = @(t) LagrangeInterp_section(t,tn,j).*exp(-lambda*(time_end-t)).*ppval(pp,t);
    a(j,1) = integral(@(t)Aj(t).*(time_end-t)/dt,time_beg,time_end);
    a(j,2) = integral(@(t)Aj(t).*(t-time_beg)/dt,time_beg,time_end);
end

zi =  lambda *( C_old*exp(-lambda*dt) )/p;
TR = TR + IV*dpdt/p;
TR = TR + lambda * ( a(end,1)*NFId_old + a(end,2)*NFId_new )/p ;

for n=1:order
    zi = zi + lambda *( a(n,1)*NFId_old*shape_beg(:,n) + a(n,2)*NFId_new*shape_beg(:,n) )/p;
end
A = IV-dt*TR;
rhs = IV*shape_beg(:,end) + dt*zi;

% apply BC
if npar.set_bc_last
    [A,rhs]=apply_BC(A,rhs,npar.add_ones_on_diagonal);
else
    rhs=apply_BC_vec_only(rhs);
end

Re = A*shape - rhs;
    
if dat.solve_prec 
    C_new =  C_old*exp(-lambda*dt) + ( a(end,1)*NFId_old + a(end,2)*NFId_new )*shape ;
    for n=1:order
        C_new = C_new + ( a(n,1)*NFId_old + a(n,2)*NFId_new )*shape_beg(:,n);
    end
    dat.C_new = C_new;
    dat.t = t;
    dat.y = y;
    dat.X = X;
end