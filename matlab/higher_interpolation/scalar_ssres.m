function dydt = scalar_ssres(time, y, a, q)

% evaluates the steady state residual. here, it is a(t).u(t)+q(t)
% I keep y as the generic solution name (y is u here).

dydt = a(time)*y + q(time);

return
end

