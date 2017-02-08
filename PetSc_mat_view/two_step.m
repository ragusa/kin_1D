close all; clear;clc

fd = PetscOpenFile('A');
A = PetscBinaryRead(fd);

fd = PetscOpenFile('D');
D = PetscBinaryRead(fd);

fd = PetscOpenFile('S');
S = PetscBinaryRead(fd);

fd = PetscOpenFile('F');
F = PetscBinaryRead(fd);

fd = PetscOpenFile('Dd');
Dd = PetscBinaryRead(fd);

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

SS = spalloc(N,N,nnz(S));
row_id=0;
for g=1:G
    for i=1:n
        ia=k(i)+(g-1);
        row=[];
        for gg=1:G
            ja=k+(gg-1);
            row = [row S(ia,ja)];
        end
        row_id=row_id+1;
        SS(row_id,:)=row;
    end
end

FF = spalloc(N,N,nnz(F));
row_id=0;
for g=1:G
    for i=1:n
        ia=k(i)+(g-1);
        row=[];
        for gg=1:G
            ja=k+(gg-1);
            row = [row F(ia,ja)];
        end
        row_id=row_id+1;
        FF(row_id,:)=row;
    end
end

DDd = spalloc(N,N,nnz(Dd));
row_id=0;
for g=1:G
    for i=1:n
        ia=k(i)+(g-1);
        row=[];
        for gg=1:G
            ja=k+(gg-1);
            row = [row Dd(ia,ja)];
        end
        row_id=row_id+1;
        DDd(row_id,:)=row;
    end
end

subplot(1,2,1)
spy(D)
subplot(1,2,2)
spy(D-Dd)