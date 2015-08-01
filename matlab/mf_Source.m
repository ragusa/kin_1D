function Source = mf_Source(imat,x,t,phi)
global dat
inv_vel = evaluate_material_prop(dat.inv_vel{imat},t,x);
nusigf_p = evaluate_material_prop(dat.nusigf_p{imat},t,x);
siga = evaluate_material_prop(dat.siga{imat},t,x);
lambda = dat.lambda;
cdiff = evaluate_material_prop(dat.cdiff{imat},t,x);
nusigf_d = evaluate_material_prop(dat.nusigf_d{imat},t,x);
L = dat.width;
if phi   
    Source = (x/L.*(1-x/L).*(inv_vel-(nusigf_p-siga+lambda)*(t+1))+2*(t+1)/L^2*cdiff);
else
    Source = x/L.*(1-x/L).*(1-(nusigf_d-lambda)*(t+1));
end