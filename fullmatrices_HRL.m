function[A,B]=fullmatrices_HRL(N,x)

% construction of matrices D^2 and D^2 - a^2 I --- x=a^2:

M=chebyder2(N,N);
X=M-x*eye(N);

% construction of bc vectors in z=0 e z=1:

v_0=chebybval(N,0);
v_1=chebybval(N,1);

% matrices initialization:

A=zeros(2*N,2*N);
B=zeros(2*N,2*N);

% construction of matrix A:

for i=1:N-2
    for j=1:N
        A(i,j)=X(i,j);
    end
end

for i=N+1:2*N-2
    for j=N+1:2*N
        A(i,j)=X(i-N,j-N);
    end
end

% complete the matrix A blocks with bc:

for j=1:N
    A(N-1,j)=v_0(j);
    A(N,j)=v_1(j);
end

for j=N+1:2*N
    A(2*N-1,j)=v_0(j-N);
    A(2*N,j)=v_1(j-N);
end

% construction of matrix B:

I=eye(N-2,N-2);

for i=1:N-2
    for j=N+1:2*N-2
        B(i,j)=-x*I(i,j-N);
    end
end

for i=N+1:2*N-2
    for j=1:N-2
        B(i,j)=-I(i-N,j);
    end
end

end