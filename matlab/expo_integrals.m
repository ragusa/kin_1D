clear all; clc
syms t1 t2 t lambda ;

A=int((t2-t)/(t2-t1)*exp(lambda*(t-t2)),t1,t2);
A=simplify(A)

B=int((t-t1)/(t2-t1)*exp(lambda*(t-t2)),t1,t2);
B=simplify(B)


clear t1 t2 lambda t
t1=1;
t2=3;
lambda=0.2;

eval(A)
eval(B)

At= @(t)(t2-t)/(t2-t1).*exp(lambda*(t-t2));
Bt= @(t)(t-t1)/(t2-t1).*exp(lambda*(t-t2));

quad(@(t)At(t),t1,t2)
quad(@(t)Bt(t),t1,t2)
