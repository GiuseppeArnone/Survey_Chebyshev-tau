function[F,err_approx_Taylor,err_approx_Cheby]=Matrix_Representation(f,N,order)
%%%% F is the symbolic function to be expanded in Chebyshev polynomials;
%%%% ORDER is the order of Taylor polynomial to approximate f in first instance;
%%%% N is the dimension of the basis of Chebyshev polynomials;

%% Taylor expansion of f

t=taylor(f,'ExpansionPoint', 1/2,'Order',order);

%% Evaluation of the error due to Taylor expansion

t=expand(t);
ff=matlabFunction(f);
tt=matlabFunction(t);
zz=0:0.01:1;
err_approx_Taylor=max(abs(ff(zz)-tt(zz)));

%% Coefficients extraction
%%%% COEFF is a vector of coefficients before powers in Taylor polynomial t

coeff=CoeffExtraction(t,order);

%% Expansion of powers of z in Chebyshev polynomials
%%%% CHEB is the matrix where entries are coefficients before Chebyshev
%%%% polynomials. r-th row refers to power r, while c-th column refers to
%%%% c-th Chebyshev polynomial
%%%% The i-th entry of T is the coefficient before the i-th Chebyshev polynomial

cheb=zeros(order,order);
T=zeros(1,order);
for i = 2:order
    for j=1:i       % recursive formula for z^(i-1)
        if (i-j)/2==fix((i-j)/2)
            T(1,j)=nchoosek(i-1,(i-j)/2);
            if j==1
                T(1,:)=T(1,:)/2;
            end
        end
    end
    T(1,:)= T(1,:)*2^(2-i);
    cheb(i,:) = T(1,:);
    T=zeros(1,order);
end
cheb(1,1)=1;

%% Costruction of VETTCOEFF
%%%% VETTCOEFF where the i-th entry is the coefficient before the
%%%% i-th Chebyshev polynomial in the expansion of f

vettcoeff=zeros(1,order);
for i = 1:order
    vettcoeff=cheb(i,:)*coeff(i)+vettcoeff;
end

%% Evaluation of the error due to Chebyshev expansion

syms x
chebpoly=chebyshevT(0:order-1,x);
f_cheb=0;
for i = 1:order
    f_cheb=chebpoly(i)*vettcoeff(i)+f_cheb; % Chebyshev expansion of f
end
gg=matlabFunction(f_cheb);
ff=matlabFunction(f);
zz=0:0.01:1;
err_approx_Cheby=max(abs(ff(zz)-gg(zz)));

%% Construction of F, matrix representation of f
vettcoeff
F=Construction_Matrix_Representation(N,vettcoeff,order);
end
