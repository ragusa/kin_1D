function T = compute_temp(Told,phi,t)

global dat

[~,n] = size(phi);
if n ~= length(t)
    error('there must be at least 2 times')
end

fissrc = zeros(size(phi));
for i=1:n
    nf = evaluate_material_prop(dat.nusigf{1},t(i),0);
    fissrc(:,i) = nf*phi(:,i)*dat.Pnorm;
end
dt = t(end)-t(1);

switch n
    case 1
        T = Told + dat.alpha*dt*fissrc;
    case 2
        T = Told + dat.alpha*dt*(fissrc(:,2) + fissrc(:,1))/2;
    otherwise
        error('Integration scheme not known')
end

return
end