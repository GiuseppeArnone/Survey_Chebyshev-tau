clear all
close all
clc

% choose the number of Chebyshev polinomials:

N=30;

% choose some TU values with associated critical wave and Rayleigh numbers:

TU=[4,6,8,10,12]

WaveC=[10.21,12.31,21.86,34.55,48.43]

Rc=[6.2080,9.9506,15.3523,21.4701,28.2444]

% select for each eigenvalue the associated eigenvector and dispose it in a matrix:

Matrix_Eig=zeros(length(TU),2*N);

% solve the generalized eigenvalue problem for fixed Wavec(k):

for k=1:length(TU)

[C,D]=fullmatrices_PPorous(N,WaveC(k),TU(k));
[EigVec,EigVal]=eig(C,D,'qz');

EigVal_pos_r=zeros(2*N,2*N);

% select real positive eigenvalues:

for i=1:2*N
    if EigVal(i,i)>0 && abs(imag(EigVal(i,i)))<10^(-10) 
        EigVal_pos_r(i,i)=EigVal(i,i);
    end
end

flag=0;

% find the index position in EigVal of the critical Rayleigh Rc(k):

for i=1:2*N
    diff=EigVal_pos_r(i,i)-Rc(k);
        if abs(diff)<10^(-4) && flag==0 
        flag=1;
        index=i;
        end
end

X=EigVal(index,index)

% copy the coefficients vector in Matrix_Eig:

Matrix_Eig(k,:)=EigVec(:,index)

end

% ________________________

% linear combine Chebyshev polynomials T_k and coefficients w_k:

syms z
Cheby_base=chebyshevT(0:N-1, z);

% W(z) in correspondence with upper temperature TU=4:

funzW1=0;
for i=1:N
    fun=Cheby_base(i)*Matrix_Eig(1,i)+funzW1;
    funzW1=fun;
end

R1 = matlabFunction(funzW1);

% W(z) in correspondence with upper temperature TU=6:

funzW2=0;
for i=1:N
    fun=Cheby_base(i)*Matrix_Eig(2,i)+funzW2;
funzW2=fun;
end

R2 = matlabFunction(funzW2);

% W(z) in correspondence with upper temperature TU=8:

funzW3=0;
for i=1:N
    fun=Cheby_base(i)*Matrix_Eig(3,i)+funzW3;
funzW3=fun;
end

R3 = matlabFunction(funzW3);

% W(z) in correspondence with upper temperature TU=10:

funzW4=0;
for i=1:N
    fun=Cheby_base(i)*Matrix_Eig(4,i)+funzW4;
funzW4=fun;
end

R4 = matlabFunction(funzW4);

% W(z) in correspondence with upper temperature TU=12:

funzW5=0;
for i=1:N
    fun=Cheby_base(i)*Matrix_Eig(5,i)+funzW5;
funzW5=fun;
end

R5 = matlabFunction(funzW5);

% ________________________

% plots of the (rescaled for TU=8,10) eigenfunctions:

x1 = 0:0.001:1;
plot(x1,R1(x1),'k','LineWidth', 1)
title('Profile of W with T_U=4')
xlabel('z');
ylabel('W(z)');

x1 = 0:0.001:1;
plot(x1,R2(x1),'k','LineWidth', 1)
title('Profile of W with T_U=6')
xlabel('z');
ylabel('W(z)');

x1 = 0:0.001:1;
plot(x1,R3(x1),'r','LineWidth', 1)
title('Non riscaled Profile of W with T_U=8')
xlabel('z');
ylabel('W(z)');

x1 = 0:0.001:1;
plot(x1,-R3(x1),'k','LineWidth', 1)
title('Profile of W with T_U=8')
xlabel('z');
ylabel('W(z)');

x1 = 0:0.001:1;
plot(x1,R4(x1),'r','LineWidth', 1)
title('Non riscaled Profile of W with T_U=10')
xlabel('z');
ylabel('W(z)');

x1 = 0:0.001:1;
plot(x1,-R4(x1),'k','LineWidth', 1)
title('Profile of W with T_U=10')
xlabel('z');
ylabel('W(z)');

x1 = 0:0.001:1;
plot(x1,R5(x1),'k','LineWidth', 1)
title('Profile of W with T_U=12')
xlabel('z');
ylabel('W(z)');

% plot of the non normalized eigenfunctions (for the sake of clarity):

x1 = 0:0.001:1;
plot(x1,R1(x1),'r',x1,R2(x1),'r',x1,-R3(x1),'r',x1,-R4(x1),'r',x1,R5(x1),'r','LineWidth', 1)
title('Unreadable Profiles of W(z)')
xlabel('z');
ylabel('W(z)');

% ________________________

% suitably normalize eigenfunctions:

% TU=4

for k=1:length(x1)
    Nw(k)=feval(R1,x1(k));
end

M1=normalize(Nw,'range');

x1 = 0:0.001:1;
plot(x1,M1,'k','LineWidth', 1)

z1=M1(1)

% TU=6

for k=1:length(x1)
    Nw(k)=feval(R2,x1(k));
end

M2=normalize(Nw,'range');

x1 = 0:0.001:1;
plot(x1,M2,'k','LineWidth', 1)

z2=M2(1)

% TU=8

for k=1:length(x1)
    Nw(k)=-feval(R3,x1(k));
end

N3=normalize(Nw,'range');

x1 = 0:0.001:1;
plot(x1,N3,'r','LineWidth', 1)

z3=N3(1)

M3=normalize(Nw,'range',[-z3 1]);

x1 = 0:0.001:1;
plot(x1,M3,'k','LineWidth', 1)

% TU=10

for k=1:length(x1)
    Nw(k)=-feval(R4,x1(k));
end

N4=normalize(Nw,'range');

x1 = 0:0.001:1;
plot(x1,N4,'r','LineWidth', 1)

z4=N4(1)

M4=normalize(Nw,'range',[-z4 1]);

x1 = 0:0.001:1;
plot(x1,M4,'k','LineWidth', 1)


% TU=12

for k=1:length(x1)
    Nw(k)=feval(R5,x1(k));
end

N5=normalize(Nw,'range');

x1 = 0:0.001:1;
plot(x1,N5,'r','LineWidth', 1)

z5=N5(1)

M5=normalize(Nw,'range',[-z5 1]);

x1 = 0:0.001:1;
plot(x1,M5,'k','LineWidth', 1)

% ________________________

% normalized eigenfunctions plot:

x1 = 0:0.001:1;
y0=zeros(1,length(x1));
plot(x1,M2,'k',x1,M3,'--k',x1,M4,'-.k',x1,M5,':k',x1,y0,'k','LineWidth', 1)
title('Profiles of W(z) normalized over the spatial layer')
xlabel('z');
ylabel('W(z)');
legend('T_U=6C','T_U=8C','T_U=10C','T_U=12C','Location',"northeast");




