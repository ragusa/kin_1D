clear;close all

fd = PetscOpenFile('A');
A = PetscBinaryRead(fd);
[N,~]=size(A);
b=ones(N,1);
restart=[];
tol=1e-8;
maxit=120;
[x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,A);
semilogy(1:length(resvec),resvec/norm(A\b));
leg = char('PJFNK+LU SMP full ltol=1e-8');
hold on;

fd = PetscOpenFile('PJFNK_GMRES_SMPfull_l8');
A1 = PetscBinaryRead(fd);
restart=[];
tol=1e-8;
maxit=120;
[x,flag,relres,iter,resvec1] = gmres(A,b,restart,tol,maxit,A1);
semilogy(1:length(resvec1),resvec1/norm(A1\b));
leg = char(leg,'PJFNK+GMRES SMP full ltol=1e-8');

fd = PetscOpenFile('PJFNK_GMRES_SMPfull_l10');
A2 = PetscBinaryRead(fd);
restart=[];
tol=1e-10;
maxit=120;
[x,flag,relres,iter,resvec2] = gmres(A,b,restart,tol,maxit,A2);
semilogy(1:length(resvec2),resvec2/norm(A2\b));
leg = char(leg,'PJFNK+GMRES SMP full ltol=1e-10');

fd = PetscOpenFile('NEWTON_GMRES_SMPfull_l8');
A3 = PetscBinaryRead(fd);
restart=[];
tol=1e-8;
maxit=120;
[x,flag,relres,iter,resvec3] = gmres(A,b,restart,tol,maxit,A3);
semilogy(1:length(resvec3),resvec3/norm(A3\b));
leg = char(leg,'NEWTON+GMRES SMP full ltol=1e-8');

fd = PetscOpenFile('PJFNK_GMRES_SMPdiag_l8');
D1 = PetscBinaryRead(fd);
restart=[];
tol=1e-8;
maxit=120;
% [x,flag,relres,iter,resvec4] = gmres(A,b,restart,tol,maxit,eye(size(D1)),D1);
[x,flag,relres,iter,resvec4] = gmres(A,b,restart,tol,maxit,D1);
semilogy(1:length(resvec4),resvec4/norm(D1\b));
leg = char(leg,'NEWTON+GMRES SMP diagonal ltol=1e-8');


legend(leg,'Location','Best')
xlabel('# of iters');
ylabel('residual norm / rhs norm ');
