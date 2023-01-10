clear all
close all
clc

% choose the number of Chebyshev polinomials and the discretization step:

N=30;
step=0.01;

% x=a^2

x=0.7:step:40; 

% choose the upper layer temperature T_U>=4

TU=8;

for i=1:length(x)

    [A,B]=fullmatrices_PPorous(N,x(i),TU);
    R=eig(A,B,'qz');

    R_pos_r=zeros(length(R),1);

% selection of real positive eigenvalues:
    
    for j=1:length(R)
        if R(j)>0 && abs(imag(R(j)))<10^(-10)
            R_pos_r(j)=R(j);
        end
    end

% ordering of positive eigenvalues:

    R_pos_r_ord=sort(R_pos_r);

% selection of the first positive eigenvalue and its square:

    flag=0;
    for j=1:length(R_pos_r_ord) 
        if R_pos_r_ord(j)>0 && flag==0
            R_x(i)=R_pos_r_ord(j);
            R_x_sqr(i)=R_x(i)^2;
            flag=1;
        end
    end

end

% determination of the critical Rayleigh number:

[M,i]=min(R_x_sqr);
WaveC=double(x(i))
Rc=R_x(i)

% select for each eigenvalue the associated eigenvector and dispose it in a matrix:

Matrix_Eig=zeros(1,2*N);

% solve the generalized eigenvalue problem for fixed Wavec and TU:

[C,D]=fullmatrices_PPorous(N,WaveC,TU);
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
    diff=EigVal_pos_r(i,i)-Rc;
        if abs(diff)<10^(-4) && flag==0 
        flag=1;
        index=i;
        end
end

X=EigVal(index,index)

% copy the coefficients vector in Matrix_Eig:

Matrix_Eig(1,:)=EigVec(:,index)

% ________________________

% linear combine Chebyshev polynomials T_k and coefficients w_k:

syms z
Cheby_base=chebyshevT(0:N-1, z);

% W(z) in correspondence with upper temperature TU:

funzW1=0;
for i=1:N
    fun=Cheby_base(i)*Matrix_Eig(1,i)+funzW1;
    funzW1=fun;
end

R = matlabFunction(funzW1);

% ________________________

% plots of the eigenfunction:

x1 = 0:0.001:1;
plot(x1,-R(x1),'k','LineWidth', 1)
title('Profile of W')
xlabel('z');
ylabel('W(z)');

% ________________________

% drow the convective roll:

wave=sqrt(WaveC)

x = 0:0.001:1;
z = 0:0.001:2*pi/wave;
[X,Z] = meshgrid(x,z);
Y = R(X).*sin(wave*Z);
figure
contour(Z,X,Y,30,'ShowText','off')
title('Convective roll perturbation')
ylabel('0<z<1');
xlabel('0<x<2pi/a_c');

% plot of w(x,z):

x = 0:0.03:1;
z = 0:0.03:2*pi/wave;
[X,Z] = meshgrid(x,z);
Y = R(X).*sin(wave*Z);
surfc(Z,X,Y)
ylabel('0<z<1');
xlabel('0<x<2pi/a_c');
zlabel('w(x,z)=W(z)sin(a_cx)');


