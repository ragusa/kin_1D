close all;

fd = PetscOpenFile('A');
A = PetscBinaryRead(fd);

fd = PetscOpenFile('D');
D = PetscBinaryRead(fd);


% % condest(A)
% % condest(D)
% % 

figure(1)
subplot(1,4,1); spy(A,5);
subplot(1,4,2); spy(D,5);

G=11;
[N,~]=size(A);
n=N/G;

k=1:G:N;
if length(k)~=n
    error('lengths');
end

AA = spalloc(N,N,nnz(A));
row_id=0;
for g=1:G
    for i=1:n
        ia=k(i)+(g-1);
        row=[];
        for gg=1:G
            ja=k+(gg-1);
            row = [row A(ia,ja)];
        end
        row_id=row_id+1;
        AA(row_id,:)=row;
    end
end
subplot(1,4,3); spy(AA,5);

DD = spalloc(N,N,nnz(D));
row_id=0;
for g=1:G
    for i=1:n
        ia=k(i)+(g-1);
        row=[];
        for gg=1:G
            ja=k+(gg-1);
            row = [row D(ia,ja)];
        end
        row_id=row_id+1;
        DD(row_id,:)=row;
    end
end
subplot(1,4,4); spy(AA-DD,5);

Sparsity = spones(AA);
sss=Sparsity(1:n,1:n);
% xs = eye(G,G);
% xs(7:G,7:G)=1;
% xs(6,7)=1;
% xs(5:8,6)=1;
% xs(2:4,1)=1;
% xs(3:4,2)=1;
% xs(4,3)=1;
% xs(5,4)=1;
% xs(6,5)=1;
% spy(kron(xs,sss)-Sparsity)

figure(99)
subplot(2,3,1); spy(A,5); title('full system')
subplot(2,3,2); spy(D,5); title('block diag precondtioner')
subplot(2,3,3); spy(A-D,5); title('full minus block diag')
subplot(2,3,4); spy(AA,5); title('full system re-ordered')
subplot(2,3,5); spy(DD,5); title('block diag precondtioner re-ordered')
subplot(2,3,6); spy(AA-DD,5); title('full minus block diag re-ordered')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


b=ones(N,1);

% matlab \ solve
t=cputime;
x_exact = A\b;
t_exact=cputime-t

% matlab default gmres w/o prec, no restart
restart=[];
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(100); semilogy(1:length(resvec),resvec/norm(b)); leg=char('gmres, no prec., no restart'); hold all;


% matlab default gmres w/o prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
semilogy(1:length(resvec),resvec/norm(b)); leg=char(leg,sprintf('gmres, no prec., restart %g',restart));

% matlab default gmres WITH prec, no restart
restart=[];
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,D);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
semilogy(1:length(resvec),resvec/norm(D\b)); leg=char(leg,'gmres, WITH prec., no restart'); 

% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,D);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
semilogy(1:length(resvec),resvec/norm(D\b)); leg=char(leg,sprintf('gmres, WITH prec., restart %g',restart));

legend(leg,'Location','Best')
title('GMRES convergence using Matlab and Rattlesnake''s matrices');
xlabel('# of iters');
ylabel('residual norm / rhs norm ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,D);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(101);
semilogy(1:length(resvec),resvec/norm(D\b)); hold all; leg=char(sprintf('gmres, WITH prec., restart %g',restart));

figure(22);v=eig(full(D\A));plot(real(v),imag(v),'o');pause;

%%%% new preconditioner
down=-1;
up=+0;
xs=tril(triu(ones(11,11),down),up);
DD=kron(xs,sss).*AA;
% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit,DD);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(101);
semilogy(1:length(resvec),resvec/norm(DD\b)); leg=char(leg,sprintf('gmres, NEW prec.(%g,%g), restart %g',down,up,restart));

figure(22);v=eig(full(DD\AA));plot(real(v),imag(v),'o');pause;

%%%% new preconditioner
down=-0;
up=+1;
xs=tril(triu(ones(11,11),down),up);
DD=kron(xs,sss).*AA;
% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit,DD);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(101);
semilogy(1:length(resvec),resvec/norm(DD\b)); leg=char(leg,sprintf('gmres, NEW prec.(%g,%g), restart %g',down,up,restart));

figure(22);v=eig(full(DD\AA));plot(real(v),imag(v),'o');pause;


%%%% new preconditioner
down=-1;
up=+1;
xs=tril(triu(ones(11,11),down),up);
DD=kron(xs,sss).*AA;
% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit,DD);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(101);
semilogy(1:length(resvec),resvec/norm(DD\b)); leg=char(leg,sprintf('gmres, NEW prec.(%g,%g), restart %g',down,up,restart));

figure(22);v=eig(full(DD\AA));plot(real(v),imag(v),'o');pause;


%%%% new preconditioner
down=-2;
up=+1;
xs=tril(triu(ones(11,11),down),up);
DD=kron(xs,sss).*AA;
% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit,DD);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(101);
semilogy(1:length(resvec),resvec/norm(DD\b)); leg=char(leg,sprintf('gmres, NEW prec.(%g,%g), restart %g',down,up,restart));

figure(22);v=eig(full(DD\AA));plot(real(v),imag(v),'o');pause;


%%%% new preconditioner
down=-1;
up=+2;
xs=tril(triu(ones(11,11),down),up);
DD=kron(xs,sss).*AA;
% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit,DD);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(101);
semilogy(1:length(resvec),resvec/norm(DD\b)); leg=char(leg,sprintf('gmres, NEW prec.(%g,%g), restart %g',down,up,restart));

figure(22);v=eig(full(DD\AA));plot(real(v),imag(v),'o');pause;


%%%% new preconditioner
down=-2;
up=+2;
xs=tril(triu(ones(11,11),down),up);
DD=kron(xs,sss).*AA;
% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit,DD);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(101);
semilogy(1:length(resvec),resvec/norm(DD\b)); leg=char(leg,sprintf('gmres, NEW prec.(%g,%g), restart %g',down,up,restart));

figure(22);v=eig(full(DD\AA));plot(real(v),imag(v),'o');pause;


%%%% new preconditioner
down=-3;
up=+3;
xs=tril(triu(ones(11,11),down),up);
DD=kron(xs,sss).*AA;
% matlab default gmres WITH prec, WITH restart
restart=15;
tol=1e-10;
maxit=120;
t=cputime;
[x,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit,DD);
t_gmres=cputime-t
if flag ~=0
    warning('gmres flag not 0: %g',flag);
    iter
    err_gmres = norm((x-x_exact)./x)
end
iter
figure(101);
semilogy(1:length(resvec),resvec/norm(DD\b)); leg=char(leg,sprintf('gmres, NEW prec.(%g,%g), restart %g',down,up,restart));

figure(22);v=eig(full(DD\AA));plot(real(v),imag(v),'o');pause;


figure(101);

legend(leg,'Location','Best')
title('GMRES convergence using Matlab and Rattlesnake''s matrices');
xlabel('# of iters');
ylabel('residual norm / rhs norm ');
