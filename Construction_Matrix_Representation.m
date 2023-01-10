function[F]=Construction_Matrix_Representation(N,vettcoeff,order)
%%%% VETTCOEFF is the vector where the i-th entry is the coefficient before
%%%% the i-th Chebyshev polynomial in the expansion of function f
%%%% ORDER is the order of Taylor polynomial to approximate f in first instance;
%%%% N is the dimension of the basis of Chebyshev polynomials;

matrix=zeros(N,order+N); % Matrix to build F
first=zeros(1,order+N); % first contribution
second=zeros(1,order+N); % second contribution
contrib=zeros(1,order+N);
for k=1:N
    for i=1:order
        first(1,i+k-1)=vettcoeff(i)/2;
        second(1,abs(i-k)+1)=vettcoeff(i)/2;
        contrib=first+second+contrib;
        first=zeros(1,order+N);
        second=zeros(1,order+N);
    end
    matrix(k,:)=contrib;
    contrib=zeros(1,order+N);
end
F=transpose(matrix);
F=F(1:N,:);
end