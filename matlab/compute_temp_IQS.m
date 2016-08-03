function T = compute_temp_IQS(Told,phi,t,p)

global dat

[~,n] = size(phi);
if length(p) ~= length(t)
    error('p and time must be same length')
elseif n~=2
    error('must have two shape functions')    
end
t1 = t(1);
t2 = t(end);
dt = t2-t1;

nf =  evaluate_material_prop(dat.nusigf{1},t1,0);
fissrc_old = nf*phi(:,1);
nf =  evaluate_material_prop(dat.nusigf{1},t2,0);
fissrc_new = nf*phi(:,2);

pp = interp1(t,p,'linear','pp');
A1 = @(t) (t2-t)/dt.*ppval(pp,t);
A2 = @(t) (t-t1)/dt.*ppval(pp,t);
a1 = integral(A1,t1,t2);
a2 = integral(A2,t1,t2);

T = Told + dat.alpha*(a2*fissrc_new + a1*fissrc_old)*dat.Pnorm;

return
end