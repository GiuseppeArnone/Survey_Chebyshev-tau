% p-q block

for i=(p-1)*N+1:p*N-2
    for j=(q-1)*N+1:q*N
        A(i,j)=X(i-(p-1)*N,j-(q-1)*N);
    end
end

% p-q bc

for j=(q-1)*N+1:q*N
    A(p*N-1,j)=bv_0(j-(q-1)*N);
    A(p*N,j)=bv_1(j-(q-1)*N);
end