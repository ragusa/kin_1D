function [S] = assemble_source_operator(curr_time)

global npar dat

n = npar.n;
S = zeros(2*n,1);

S_phi = assemble_source(dat.source_phi,curr_time);
S_c   = assemble_source(dat.source_c,curr_time);

S(1:n)     = S_phi;
S(n+1:end) = S_c;

return
end