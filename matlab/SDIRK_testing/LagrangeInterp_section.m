function P = LagrangeInterp_section(t,tn,j)

P = 1.0;
for k=1:length(tn)
    if k ~= j
        P = P .* (t-tn(k))/(tn(j)-tn(k));
    end
end