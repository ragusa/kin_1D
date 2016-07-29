function T = compute_temp(Told,phi,t)

global dat

[~,n] = size(phi);
if n ~= length(t)
    error('phi and time must be same length')
end

fissrc = zeros(size(phi));
for i=1:n
    NF = assemble_mass(dat.nusigf,t(i));
    fissrc(:,i) = NF*phi(:,i)*dat.Pnorm;
end
dt = t(2)-t(1);

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