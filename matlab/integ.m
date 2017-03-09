syms t1 t2 L t

int(exp(L*t), t1, t2)


A=int((t2-t)/(t2-t1)*exp(L*t), t1, t2)
B=expand(exp(-L*t2)*A);
simplify(B, 'IgnoreAnalyticConstraints', true)
% (exp(-L*t2)*(exp(L*t1) - exp(L*t2) - L*t1*exp(L*t1) + L*t2*exp(L*t1)))/(L^2*(t1 - t2))

C=int((t-t1)/(t2-t1)*exp(L*t), t1, t2)
D=expand(exp(-L*t2)*C);
simplify(D, 'IgnoreAnalyticConstraints', true)
% -(exp(L*(t1 - t2)) - L*t1 + L*t2 - 1)/(L^2*(t1 - t2))

simplify(B+D)
% -(exp(L*(t1 - t2)) - 1)/L