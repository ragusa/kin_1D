function S = assemble_Source_vector(curr_time)

global npar

n   = npar.n;
S   = zeros(2*n, 1);

% flux-flux
S_phi = assemble_source(@mf_Source,curr_time,true);
S_C = assemble_source(@mf_Source,curr_time,false);

S(1:n) = S_phi;
S(n+1:2*n) = S_C;

return
end