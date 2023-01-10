clear all
close all
clc

N=30;

WaveC=4.9348
Rc=25.6420

% select for each eigenvalue the associated eigenvector and dispose it in a matrix:

Matrix_Eig=zeros(1,3*N);

% solve the generalized eigenvalue problem for fixed Wavec and TU:

[C,D]=fullmatrices_B(N,WaveC);
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
plot(x1,R(x1),'k','LineWidth', 1)
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
sc=surfc(Z,X,Y)
ylabel('0<z<1');
xlabel('0<x<2pi/a_c');
zlabel('w(x,z)=W(z)sin(a_cx)');
